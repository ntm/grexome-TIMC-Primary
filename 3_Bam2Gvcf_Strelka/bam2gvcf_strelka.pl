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



#############################################
## options / params from the command-line

my $USAGE = '
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--outdir string [required] : dir where GVCF-containing subdirs will be
    created (one subdir per grexome)
--indir string [defaults to a working path on krakenator, luxor and fauve] : 
    subdir containing the BAMs
--strelka [default "strelkaGermline.sh"] : must be either the path+name of
    configureStrelkaGermlineWorkflow.py (from strelka distrib), or a wrapper script
    (eg strelkaGermline.sh on fauve, this is the default) that can be in your PATH
--first int [default = 50] : first grexomeNum to process (>= 50)
--last int [default = 9999] : last grexomeNum to process (>= $first)
--genome string [defaults to searching in NTM\'s standard places] : path+filename
    of the reference genome in fasta format, subdir must also contain chroms-only BED
---jobs N [default = 4] : number of cores that strelka can use
--real : actually do the work, otherwise this is a dry run, just print 
    info on what would be done
--help : print this USAGE';


# dir where GVCF-containing subdirs will be created, no default => MUST BE PROVIDED
my $outDir = "";

# subdir where BAMs are stored, defaults to current correct dir on 
# luxor (where the BAMS are currently stored), and also on fauve (and
# others) via ~/sshMounts/luxor/
my $inDir = "PierreRay_DATA/BAMs_All_Selected/";

if (-e "/home/nthierry/sshMounts/luxor/") {
    $inDir = "/home/nthierry/sshMounts/luxor/$inDir";
}
else {
    $inDir = "/home/nthierry/$inDir";
}

# first and last grexomeNums to process, default to everything!
my ($firstGrex, $lastGrex) = (50,9999);

# path+name of configureStrelkaGermlineWorkflow.py from strelka distrib,
# or use a one-liner wrapper eg strelkaGermline.sh that can be in your PATH,
# this is the default (works on fauve, luxor...)
my $strelka = "strelkaGermline.sh";

# full path to ref genome, must be the one used for producing the BAMs
my $refGenome;

# number of parallel jobs to run, get from command-line --jobs, default to 4
my $jobs = 4;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

GetOptions ("outdir=s" => \$outDir,
	    "indir=s" => \$inDir,
	    "first=i" => \$firstGrex,
	    "last=i" => \$lastGrex,
	    "strelka=s" => \$strelka,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) &&
    die "$USAGE\n\n";
($outDir) || 
    die "$USAGE\n\nE: you MUST specify the dir where GVCF-containing subdirs will be created, with --outdir\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E: outDir $outDir doesn't exist as a dir but can't be created\n";
(-d $inDir) || 
    die "E: inDir $inDir doesn't exist or isn't a dir\n";
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

# make sure refGenome exists if provided
if ($refGenome) {
    (-f "$refGenome") ||
	die "E: the provided reference genome $refGenome doesn't exist or can't be read\n";
}
else {
    # if refGenome wasn't provided try to find it in several default subdirs
    # so it works on fauve, luxor, luke, dahu...
    if (-f "/data/HumanGenome/hs38DH.fa") {
	$refGenome = "/data/HumanGenome/hs38DH.fa";
    }
    elsif (-f "/home/nthierry/HumanGenome/hs38DH.fa") {
	$refGenome = "/home/nthierry/HumanGenome/hs38DH.fa";
    }
    elsif (-f "/bettik/nthierry/HumanGenome/hs38DH.fa") {
	$refGenome = "/bettik/nthierry/HumanGenome/hs38DH.fa";
    }
    else {
	die "E: cannot find human ref genome in the default paths, you must use eg --refgenome path/to/ref/hs38DH.fa";
    }
}

# full path to hs38_chrom-only BED file, to ignore decoy/unplaced/alt 
# regions, should be alongside the ref genome
my $chromsBed = $refGenome;
($chromsBed =~ s/hs38DH.fa$/hs38_chroms.bed.gz/) ||
    die "E: cannot substitute hs38 fasta for bed\n";
(-f $chromsBed) ||
    die "E: chromsBed file $chromsBed doesn't exist, rsync it from somewhere or fix the code";



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
