#!/usr/bin/perl


# 16/06/2019
# NTM

# Call variants using strelka on BAM files.

# GVCF files and strelka stats will be created in subdirs of $outDir,
# if a subdir matching an inFile already exists this inFile is skipped.
# When anything potentially important happens (eg skipping an infile),
# messages are logged on stdout.
#
# See $USAGE for arguments.

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename qw(basename);
use FindBin qw($RealBin);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## options / params from the command-line


# subdir where BAMs can be found
my $inDir;

# dir where GVCF-containing subdirs will be created
my $outDir;

# first and last grexomeNums to process, default to everything!
my ($firstGrex, $lastGrex) = (50,9999);

# path+name of configureStrelkaGermlineWorkflow.py from strelka distrib,
# or use a one-liner wrapper eg strelkaGermline.sh that can be in your PATH,
# this is the default (works on fauve, luxor...)
my $strelka = "strelkaGermline.sh";

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/../grexomeTIMCprim_config.pm";

# number of parallel jobs to run, get from command-line --jobs, default to 4
my $jobs = 12;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir string [no default] : subdir containing the BAMs
--outdir string [no default] : dir where GVCF-containing subdirs will be
    created (one subdir per grexome)
--first int [$firstGrex] : first grexomeNum to process (>= 50)
--last int [$lastGrex] : last grexomeNum to process (>= first)
--strelka [default \"$strelka\"] : must be either the path+name of
    configureStrelkaGermlineWorkflow.py (from strelka distrib), or a wrapper script
    (eg strelkaGermline.sh) that can be in your PATH
---jobs N [default = $jobs] : number of cores that strelka can use
--real : actually do the work, otherwise this is a dry run, just print 
    info on what would be done
--help : print this USAGE";


GetOptions ("indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "first=i" => \$firstGrex,
	    "last=i" => \$lastGrex,
	    "strelka=s" => \$strelka,
	    "config=s" => \$config,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n$USAGE\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenome refGenomeChromsBed) );

($inDir) ||
    die "E $0: you MUST provide a dir where BAMs can be found, with --indir\n$USAGE\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

($outDir) || 
    die "$USAGE\n\nE: you MUST specify the dir where GVCF-containing subdirs will be created, with --outdir\n$USAGE\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E: outDir $outDir doesn't exist as a dir but can't be created\n";

(($firstGrex >= 50) && ($lastGrex >= $firstGrex)) ||
    die "E: first grexomeNum must be >=50 and last must be >=first\n";

# bring $lastGrex down to the largest existing grexome*.bam in inDir
# (mostly in case we have the default 9999)
while($lastGrex > $firstGrex) {
    my $grexome = $lastGrex;
    # left-pad with zeroes to 4 digits
    ($grexome < 10) && ($grexome = "0$grexome");
    ($grexome < 100) && ($grexome = "0$grexome");
    ($grexome < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";
    (-f "$inDir/${grexome}.bam") && last;
    $lastGrex--;
}

# make sure strelka can be found
(`which $strelka` =~ /$strelka$/) || 
    die "E: the strelka python (or shell wrapper) $strelka can't be found\n";

# ref genome and BED with chromosomes 1-22, X, Y, M
my $refGenome = &refGenome();
my $chromsBed = &refGenomeChromsBed();


#############################################
## process each sample of interest

my $now = strftime("%F %T", localtime);
print "I: $now - $0 STARTING TO WORK\n";

foreach my $gNum ($firstGrex..$lastGrex) {
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";

    # make sure we have bam and bai files for $grexome, otherwise skip
    my $bam = "$inDir/${grexome}.bam";
    ((-e $bam) && ((-e "$bam.bai") || (-e "$inDir/${grexome}.bai"))) ||
	((print "W: $grexome was asked for but we don't have a BAM and BAI for it in $inDir, skipping this grexome\n") && next);

    # strelka will produce a GVCF but also other files (stats etc) in $runDir
    my $runDir = "$outDir/$grexome/";

    # strelka configure command
    my $com = "$strelka --bam $bam  --referenceFasta $refGenome --callRegions $chromsBed --exome --runDir $runDir";

    # only run the strelka configuration step if runDir doesn't exist
    if (! -e $runDir) {
	if (! $real) {
	    print "I: dryrun, would configure strelka for $grexome with: $com\n";
	}
	else {
	    $now = strftime("%F %T", localtime);
	    print "I: $now - configuring strelka for $grexome with command: $com\n";
	    system($com);
	}
    }
    else {
	print "I: runDir $runDir already exists, assuming strelka is already configured\n";
    }

    # now run strelka (does nothing if it was already completed, but resumes 
    # if it was interrupted, nice)
    # NOTE: we use --quiet since strelka is very verbose and stderr is replicated
    # to workspace/pyflow.data/logs/pyflow_log.txt anyways
    $com = "$runDir/runWorkflow.py -m local -j $jobs --quiet";
    if (! $real) {
	print "I: dryrun, would run strelka for $grexome with: $com\n";
    }
    else {
	(-e "$runDir/runWorkflow.py") ||
	    ((print "I: want to run strelka for $grexome but runDir $runDir doesn't exist, configuring probably failed\n") && next);
	$now = strftime("%F %T", localtime);
	print "I: $now - running strelka for $grexome with command: $com\n";
	system($com);
    }
}

$now = strftime("%F %T", localtime);
print "I: $now - $0 ALL DONE\n";
