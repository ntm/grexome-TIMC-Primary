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

# comma-separated list of samples (FASTQs) to process (required)
my $samples = '';

# dir where GVCF-containing subdirs will be created
my $outDir;

# path+name of configureStrelkaGermlineWorkflow.py from strelka distrib,
# or use a one-liner wrapper eg strelkaGermline.sh that can be in your PATH,
# this is the default (works on fauve, luxor...)
my $strelka = "strelkaGermline.sh";

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCprim_config.pm";

# number of parallel jobs to run
my $jobs = 16;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir : subdir containing the BAMs
--samples : comma-separated list of sampleIDs to process, for each sample we expect
	  [sample].bam and [sample].bam.bai files in indir
--outdir : dir where GVCF-containing subdirs will be created (one subdir per sample)
--strelka [$strelka] : must be either the path+name of configureStrelkaGermlineWorkflow.py
          (from strelka distrib), or a wrapper script that can be in your PATH
--config [$config] : your customized copy (with path) of the distributed *config.pm
--jobs [$jobs] : number of cores that strelka can use
--real : actually do the work, otherwise this is a dry run
--help : print this USAGE";


GetOptions ("indir=s" => \$inDir,
	    "samples=s" => \$samples,
	    "outdir=s" => \$outDir,
	    "strelka=s" => \$strelka,
	    "config=s" => \$config,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$0 $USAGE\n\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n$USAGE\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenome refGenomeChromsBed) );

($inDir) ||
    die "E $0: you MUST provide --indir where BAMs can be found\n$USAGE\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

# save samples in %samples to detect duplicates and allow sorting
my %samples;
foreach my $sample (split(/,/, $samples)) {
    if ($samples{$sample}) {
	warn "W $0: sample $sample was specified twice, is that a typo? Ignoring the dupe\n";
	next;
    }
    $samples{$sample} = 1;
}

($outDir) || 
    die "E $0: you MUST specify --outdir where GVCF-containing subdirs will be created\n$USAGE\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E $0: outDir $outDir doesn't exist as a dir and can't be created\n";

# make sure strelka can be found
(`which $strelka` =~ /$strelka$/) || 
    die "E $0: the strelka python (or shell wrapper) $strelka can't be found\n";

#############################################

# ref genome and BED with chromosomes 1-22, X, Y, M
my $refGenome = &refGenome();
my $chromsBed = &refGenomeChromsBed();

#############################################
## process each sample of interest

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - STARTING TO WORK\n";

foreach my $sample (sort keys(%samples)) {
    # make sure we have bam and bai files for $sample, otherwise skip
    my $bam = "$inDir/$sample.bam";
    ((-e $bam) && (-e "$bam.bai")) || 
	((warn "W $0: no BAM or BAI for $sample in inDir $inDir, skipping $sample\n") && next);

    # strelka will produce a GVCF but also other files (stats etc) in $runDir
    my $runDir = "$outDir/$sample/";

    # strelka configure command
    my $com = "$strelka --bam $bam  --referenceFasta $refGenome --callRegions $chromsBed --exome --runDir $runDir > /dev/null";

    # only run the strelka configuration step if runDir doesn't exist
    if (! -e $runDir) {
	if (! $real) {
	    warn "I $0: dryrun, would configure strelka for $sample with: $com\n";
	}
	else {
	    $now = strftime("%F %T", localtime);
	    warn "I $now: $0 - configuring strelka for $sample\n";
	    system($com);
	}
    }
    else {
	warn "I $0: runDir $runDir already exists, assuming strelka is already configured\n";
    }

    # now run strelka (does nothing if it was already completed, but resumes 
    # if it was interrupted, nice)
    # NOTE: we use --quiet since strelka is very verbose and stderr is replicated
    # to workspace/pyflow.data/logs/pyflow_log.txt anyways
    $com = "$runDir/runWorkflow.py -m local -j $jobs --quiet";
    if (! $real) {
	warn "I $0: dryrun, would run strelka for $sample with: $com\n";
    }
    else {
	(-e "$runDir/runWorkflow.py") ||
	    ((warn "I $0: want to run strelka for $sample but runDir $runDir doesn't exist, configure failed??\n") && next);
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - running strelka for $sample\n";
	system($com);
    }
}

$now = strftime("%F %T", localtime);
warn "I $now: $0 - ALL DONE\n";
