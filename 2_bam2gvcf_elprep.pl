#!/usr/bin/perl


# 05/07/2021
# NTM


# Call variants on BAMs and produce GVCF files using elPrep5.
#
# See $USAGE for arguments.

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename qw(basename);
use File::Temp qw(tempdir);
use FindBin qw($RealBin);
use Parallel::ForkManager;


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# elprep offers 2 possible modes:
# - filter needs a LOT of RAM but should run faster
# - sfm == split-filter-merge splits by chromosome, less RAM, more temp HDD, slower
# In filter mode we process samples sequentially (to avoid the OOM-killer),
# while in sfm mode we run up to min($jobs, $#samples) elprep jobs in parallel.
my $mode = "sfm";

# subdir where BAMs can be found
my $inDir;

# comma-separated list of samples (FASTQs) to process (required)
my $samples = '';

# dir where GVCFs will be created
my $outDir;

# dir for logs
my $logDir = "logs_elPrep/";

# path+name of elPrep5 binary, defaults to "elprep" which should be
# in your PATH
my $elprep = "elprep";

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCprim_config.pm";

# number of available cores:
# - in sfm mode we process up to $jobs samples in parallel and each process
#   gets as many threads as possible to aim for ~ $jobs threads in total;
# - in filter mode we run elprep jobs sequentially, each with $jobs threads.
my $jobs = 16;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--mode [$mode] : elprep mode, filter or sfm (ie split-filter-merge)
--indir : subdir containing the BAMs
--samples : comma-separated list of sampleIDs to process, for each sample we expect
	  [sample].bam and [sample].bam.bai files in indir
--outdir : dir where GVCF files will be created
--logdir [$logDir] : dir where elprep logs will be created
--elprep [default to \"$elprep\" which should be in PATH] : full path to elprep executable
--config [$config] : your customized copy (with path) of the distributed *config.pm
--jobs [$jobs] : number of cores/threads/jobs that we can use
--real : actually do the work, otherwise this is a dry run
--help : print this USAGE";


GetOptions ("mode=s" => \$mode,
	    "indir=s" => \$inDir,
	    "samples=s" => \$samples,
	    "outdir=s" => \$outDir,
	    "logdir=s" => \$logDir,
	    "elprep=s" => \$elprep,
	    "config=s" => \$config,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n$USAGE\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenomeElPrep refGenomeChromsBed fastTmpPath) );

($mode eq 'filter') || ($mode eq 'sfm') ||
    die "E $0: illegal --mode, you must use filter or sfm\n$USAGE\n";

($inDir) ||
    die "E $0: you MUST provide --indir where BAMs can be found\n$USAGE\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

# save samples in %samples to detect duplicates and allow sorting
my %samples;
my $nbSamples = 0;
foreach my $sample (split(/,/, $samples)) {
    if ($samples{$sample}) {
	warn "W $0: sample $sample was specified twice, is that a typo? Ignoring the dupe\n";
	next;
    }
    $samples{$sample} = 1;
    $nbSamples++;
}

($outDir) || 
    die "E $0: you MUST specify --outdir where GVCFs will be created\n$USAGE\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E $0: outDir $outDir doesn't exist as a dir and can't be created\n";

# make sure elprep executable is found
(`which $elprep` =~ /$elprep$/) ||
    die "E $0: cannot find elprep binary, you must provide it with --elprep\n";


#############################################

# ref genome
my $refGenome = &refGenomeElPrep();
# BED with chromosomes 1-22, X, Y, M
my $chromsBed = &refGenomeChromsBed();

# elprep threads for each sample, and number of samples to process in parallel
my ($threadsPerSample, $samplesInPara) = ($jobs,1);
if ($mode eq 'sfm') {
    $threadsPerSample = int (0.5 + $jobs / $nbSamples);
    # at least one thread, even if tons of samples
    ($threadsPerSample) || ($threadsPerSample=1);
    $samplesInPara = $jobs;
    ($nbSamples < $jobs) && ($samplesInPara = $nbSamples);
}

#############################################
## build generic elPrep command-line:
# -> $cmdStart $bamIn $bamOut --haplotypecaller $gvcf $cmdEnd

my $cmdStart = "$elprep $mode ";
my $cmdEnd = " --nr-of-threads $threadsPerSample ";
$cmdEnd .= "--reference $refGenome ";
$cmdEnd .= "--log-path $logDir ";
# limit to regular chromosomes (no decoy, unmapped, HLA etc)
$cmdEnd .= "--target-regions $chromsBed ";
# $cmdEnd .= "--reference-confidence GVCF " ## this is the default
# $cmdEnd .= "--timed "; ## not really useful for us

# tmpDir: only in sfm mode, elprep dies if we use --tmp-path in filter mode
if ($mode eq 'sfm') {
    my $tmpDir = tempdir(DIR => &fastTmpPath(), CLEANUP => 1);
    $cmdEnd .= "--tmp-path $tmpDir ";
}

# GQ bands (for GVCF REF blocks): too many by default, we
# want larger bands but can't be specified yet with elprep, see
# https://github.com/ExaScience/elprep/issues/52
# { 
#     my $gqb = 5;
#     while ($gqb < 80) {
# 	$cmdEnd .= " -GQB $gqb";
# 	if ($gqb < 30) { $gqb += 3;}
# 	elsif ($gqb < 50) { $gqb += 5;}
# 	else {$gqb += 10;}
#     }
# }

#############################################
## call variants

my $pm = new Parallel::ForkManager($samplesInPara);
{
    my $now = strftime("%F %T", localtime);
    warn "I $now: $0 - STARTING TO WORK, mode $mode processing $samplesInPara samples in parallel each with $threadsPerSample threads\n";
}

foreach my $sample (sort keys(%samples)) {
    # make sure we have bam and bai files for $sample, otherwise skip
    my $bam = "$inDir/$sample.bam";
    ((-e $bam) && (-e "$bam.bai")) || 
	((warn "W $0: no BAM or BAI for $sample in inDir $inDir, skipping $sample\n") && next);

    # gvcf to produce
    my $gvcf = "$outDir/${sample}.g.vcf.gz";
    # don't squash existing outfiles
    (-e "$gvcf") && 
	(warn "W $0: GVCF for $sample already exists in outDir $outDir, skipping $sample\n") && next;

    # OK build the full command, sending bamOut to /dev/null and also
    # discarding stderr (it gets duplicated in $logDir)
    my $fullCmd = "$cmdStart $bam /dev/null --haplotypecaller $gvcf $cmdEnd 2> /dev/null";

    if (! $real) {
        warn "I $0: dryrun, would run elPrep5-HaplotypeCaller for $sample with:\n$fullCmd\n";
    }
    else {
	$pm->start && next;
	my $now = strftime("%F %T", localtime);
	warn "I $now: $0 - starting elPrep5-HaplotypeCaller for $sample\n";
        if (system($fullCmd) != 0) {
            $now = strftime("%F %T", localtime);
            warn "E $now: $0 - running elPrep5-HaplotypeCaller for $sample FAILED ($?)! INSPECT THE LOGS IN $logDir\n";
        }
	else{
	    $now = strftime("%F %T", localtime);
	    warn "I $now: $0 - running elPrep5-HaplotypeCaller for $sample completed successfully\n";
	}
        $pm->finish;
    }
}
$pm->wait_all_children;

{
    my $now = strftime("%F %T", localtime);
    warn "I $now: $0 - ALL DONE\n";
}

