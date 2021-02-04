#!/usr/bin/perl

# 13/06/2020
# Initial version by OB, largely inspired by bam2gvcf_strelka.pl
# Further work by NTM


# Call variants on BAMs and produce GVCF files using GATK4.
# Also produce one log-file per sample with GATK's log,
# this script prints its own logs to stderr.


use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Parallel::ForkManager;


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

# path+name of GATK wrapper distributed with GATK4, defaults to "gatk"
# which should be in your PATH
my $gatk = "gatk";

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/../grexomeTIMCprim_config.pm";

# number of parallel jobs to run.
# NOTE: GATK4 is single-threaded and it doesn't seem possible to run it
# multi-threaded! So we run it in parallel on $jobs samples.
# Omar benchmarked this on luxor: in his hands, 
# $jobs==1 => ~2hours/sample, $jobs==15 => ~20min/sample.
my $jobs = 15;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir string [no default] : subdir containing the BAMs
--outdir string [no default] : dir where GVCF files will be created
--first int [$firstGrex] : first grexomeNum to process (>= 50)
--last int [$lastGrex] : last grexomeNum to process (>= first)
--gatk [default to \"$gatk\" which should be in PATH] : full path to gatk executable
--config string [$config] : your customized copy (with path) of the distributed *config.pm
--jobs N [default = $jobs] : number of samples to process in parallel.
--real : actually do the work, otherwise this is a dry run, just print 
    info on what would be done
--help : print this USAGE";


GetOptions ("indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "first=i" => \$firstGrex,
	    "last=i" => \$lastGrex,
	    "gatk=s" => \$gatk,
	    "config=s" => \$config,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n$USAGE\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenome refGenomeChromsBed fastTmpPath) );

($inDir) ||
    die "E $0: you MUST provide a dir where BAMs can be found, with --indir\n$USAGE\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

($outDir) || 
    die "E $0: you MUST specify the dir where GVCF-containing subdirs will be created, with --outdir\n$USAGE\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E $0: outDir $outDir doesn't exist as a dir but can't be created\n";

(($firstGrex >= 50) && ($lastGrex >= $firstGrex)) ||
    die "E $0: first grexomeNum must be >=50 and last must be >=first\n";

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

# make sure gatk executable is found, this test is disabled if
# we will be running GATK from a singularity container
($gatk =~ /singularity/) ||
    (`which $gatk` =~ /$gatk$/) ||
    die "E $0: cannot find 'gatk' (from GATK4 package), you must provide it with --gatk, you provided:\n$gatk\n";

# ref genome and BED with chromosomes 1-22, X, Y, M
my $refGenome = &refGenome();
my $chromsBed = &refGenomeChromsBed();

# tmp dir
my $tmpDir = &fastTmpPath();


my $pm = new Parallel::ForkManager($jobs);
{
    my $now = strftime("%F %T", localtime);
    warn "I: $now - $0 STARTING TO WORK\n";
}

#############################################
## build the generic GATK command-line common for all samples

# TODO: compare runtimes with larger/smaller -Xmx ?
my $cmd = "$gatk --java-options \"-Xmx8g\" HaplotypeCaller";
$cmd .= " -R $refGenome --emit-ref-confidence GVCF";

# GQ bands (for GVCF REF blocks): too many by default, we
# want slightly larger bands
{ 
    my $gqb = 5;
    while ($gqb < 80) {
	$cmd .= " -GQB $gqb";
	if ($gqb < 30) { $gqb += 3;}
	elsif ($gqb < 50) { $gqb += 5;}
	else {$gqb += 10;}
    }
}

# default loglevel is INFO but GATK is very noisy. However a lot of noise
# is WARN-level, while some interesting stuff is INFO-level...
# so for now we keep the default INFO level
# $cmd .= " --verbosity WARNING";

# in case we stay in INFO mode: log progress every 10m not 10s,
# a run can take 4h!
$cmd .= " --seconds-between-progress-updates 600";

# Omar made some tests with --native-pair-hmm-threads but 
# concluded the default 4 was ok (he saw no speedup with larger values)

# limit to regular chromosomes (no decoy, unmapped, HLA etc)
$cmd .= " -L $chromsBed";

# don't know why GATK doesn't default to FASTEST for SW
$cmd .= " --smith-waterman FASTEST_AVAILABLE";

# use fast temp storage
$cmd .= " --tmp-dir $tmpDir";

# -G -A -AX : annotations to add or exclude, defaults for now


#############################################
## call variants

foreach (my $gNum = $firstGrex; $gNum<=$lastGrex; $gNum++) {
    # grexome name
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";

    # make sure we have bam and bai files for $grexome, otherwise skip
    my $bam = "$inDir/${grexome}.bam";
    ((-e $bam) && (-e "$bam.bai")) || 
	((warn "W $0: no BAM or BAI for $grexome in inDir $inDir, skipping $grexome\n") && next);

    # gvcf to produce
    my $gvcf = "$outDir/${grexome}.g.vcf.gz";
    # don't squash existing outfiles
    (-e "$outDir/${grexome}.g.vcf.gz") && 
	(warn "W $0: GVCF for $grexome already exists in outDir $outDir, skipping $grexome.\n") && next;
    
    # OK build the full GATK command
    my $fullCmd = "$cmd -I $bam -O $gvcf";
    # GATK logging: one file per sample
    my $log = "$outDir/${grexome}.log";
    $fullCmd .= " &> $log";
    # if running via singularity assume --gatk contained opening '"(' and 
    # close them here (this is ugly but I can't figure out a cleaner way
    # to capture stderr when run via singularity)
    ($gatk =~ /singularity/) && ($fullCmd .= ' ) " ');

    if (! $real) {
        warn "I $0: dryrun, would run GATK4-HaplotypeCaller for $grexome with:\n$fullCmd\n";
    }
    else {
	$pm->start && next;
	my $now = strftime("%F %T", localtime);
	warn "I: $now - $0 starting GATK4-HaplotypeCaller for $grexome\n";
        if (system($fullCmd) != 0) {
            $now = strftime("%F %T", localtime);
            warn "E: $now - $0 running GATK4-HaplotypeCaller for $grexome FAILED ($?)! INSPECT THE LOGFILE $log\n";
        }
	else{
	    $now = strftime("%F %T", localtime);
	    warn "I: $now - $0 running GATK4-HaplotypeCaller for $grexome completed successfully\n";
	}
        $pm->finish;
    }
}
$pm->wait_all_children;

{
    my $now = strftime("%F %T", localtime);
    warn "I: $now - $0 ALL DONE\n";
}
