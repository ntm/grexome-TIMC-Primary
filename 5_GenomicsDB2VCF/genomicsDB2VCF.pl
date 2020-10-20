#!/usr/bin/perl

# 17/10/2020
# NTM


# Perform so-called "joint genotyping" with GATK GenotypeGVCFs
# from a GenomicsDB.


use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use POSIX qw(strftime);
use File::Basename qw(basename);
use FindBin qw($RealBin);
use File::Temp qw(tempdir);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## options / params from the command-line

# subdir holding the GenomicsDB
my $inDir;

# output file to create, with path.
# path must exist, and we actually produce a bgzipped VCF,
# adding .vcf.gz or .gz to $outFile if needed.
# We also create subdir [$outFile-with-no-extension]_GATKlogs/ for the GATK logs.
my $outFile;

# path+name of GATK wrapper distributed with GATK4, defaults to "gatk"
# which should be in your PATH
my $gatk = "gatk";

# bcftools binary, with path if needed, default to looking in PATH
my $bcftools = "bcftools";

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/../grexomeTIMCprim_config.pm";

# number of parallel jobs to run. GATK4 GenotypeGVCFs is single-threaded
# and very slow, we parallelize by hand running separately on each chrom
# and then merging the resulting VCFs.
my $jobs = 12;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = 'Produce VCF from a GenomicsDB by joint genotyping.
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir string [no default] : subdir holding a GenomicsDB
--outfile string [no default] : bgzipped vcf file to produce, with path (path must exist, 
     .gz or .vcf.gz are added if absent, and logdir is created based on outfile)
--gatk [default to "gatk" which could be in PATH] : full path to gatk executable
--bcftools [default to "bcftools" which could be in your PATH] : bcftools executable (for merging)
--config string [$config] : your customized copy (with path) of the distributed *config.pm
--jobs N [default = $jobs] : max number of parallel jobs to run
--real : actually do the work, otherwise this is a dry run, just print info on what would be done
--help : print this USAGE';


GetOptions ("indir=s" => \$inDir,
	    "outfile=s" => \$outFile,
	    "gatk=s" => \$gatk,
	    "bcftools=s" => \$bcftools,
	    "config=s" => \$config,
	    "jobs=i" => \$jobs,
	    "real" => \$real,
	    "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenome refGenomeChromsBed fastTmpPath) );

($inDir) ||
    die "E $0: you MUST provide a dir holding a GenomicsDB, with --indir\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

($outFile) || 
    die "E $0: you MUST specify path/to/VCF to produce, with --outfile\n";
# strip .gz and .vcf extensions if present
$outFile =~ s/\.gz$//;
$outFile =~ s/\.vcf$//;
# GATK logs will go in $logDir, squashing any pre-existing files
my $logDir = $outFile."_GATKlogs/";
(-d $logDir) || (mkdir($logDir)) || 
    die "E $0: logDir $logDir doesn't exist and cannot be mkdir'd\n";

# make sure gatk executable is found, this test is disabled if
# we will be running GATK from a singularity container
($gatk =~ /singularity/) ||
    (`which $gatk` =~ /$gatk$/) ||
    die "E $0: cannot find 'gatk' (from GATK4 package), you must provide it with --gatk, you provided:\n$gatk\n";

# make sure bcftools binary is found (crude test, might have to make it smarter)
(`which $bcftools` =~ /bcftools$/) ||
    die "E $0: cannot find bcftools, you must provide path/to/bcftools with --bcftools and program must be named 'bcftools'";


# Create subdir in &fastTmpPath so we can CLEANUP when we
# are done, because GATK leaves a lot of sh*t behind.
my $tmpDir = tempdir(DIR => &fastTmpPath(), CLEANUP => 1);

# Also make a temp subdir of $logDir where we will create the
# temporary chrom-specific VCFs (not placing them in $tmpDir
# because these VCFs are huge while &fastTmpPath can be on a 
# smallish ramdisk)
my $tmpVcfDir = tempdir(DIR => $logDir, CLEANUP => 1);

# ref genome
my $genome = &refGenome();

# construct array of chroms so we can parallelize
my @chroms;
# also construct array of corresponding temporary VCF filenames
my @tempVcfs;
# grab first column of $chromsBed == gzipped BED with chromosomes 1-22, X, Y, M
my $chromsBed = &refGenomeChromsBed();
open(CHROMS, "gunzip -c $chromsBed |") ||
    die "E $0: cannot gunzip-open chromsBed $chromsBed : $!\n";
while(my $line = <CHROMS>) {
    chomp($line);
    ($line =~ /^([^\t]+)\t/) ||
	(warn "W $0: cannot extract chrom from chromsBed $chromsBed line, skipping this line:\n$line\n" && next);
    my $chr = $1;
    push(@chroms,$chr);
    my $out = "$tmpVcfDir/tmpvcf_$chr.vcf";
    push(@tempVcfs,$out);
}
close(CHROMS);


my $pm = new Parallel::ForkManager($jobs);
{
    my $now = strftime("%F %T", localtime);
    warn "I: $now - $0 STARTING TO WORK, CALLING VARIANTS FROM $inDir\n";
}

#############################################
# build the generic GATK command-line for all chroms

my $cmd = "$gatk --java-options \"-Xmx8g\" GenotypeGVCFs";
$cmd .= " -R $genome -V gendb://$inDir";

# GATK is very noisy, quiet it down
$cmd .= " --seconds-between-progress-updates 600";

# don't create vcf.idx, if anything we will bgzip and tabix-index the VCF
$cmd .= " --create-output-variant-index false";

# use fast tmp storage
$cmd .= " --tmp-dir $tmpDir";

# -G -A -AX : annotations to add or exclude, defaults for now
# --genomicsdb-use-bcf-codec (should be faster but has some unclear limitations)


#############################################
# call variants in parallel for each chrom
foreach my $chri (0..$#chroms) {
    # limit to this chrom
    my $fullCmd = "$cmd --intervals $chroms[$chri]";
    # output to temp VCF for this chrom
    $fullCmd .= " -O $tempVcfs[$chri]";
    # chrom-specific logfile
    my $log = "$logDir/$chroms[$chri].log";
    $fullCmd .= " &> $log";


    if (! $real){
	warn "I $0: dryrun, would run GATK4-GenotypeGVCFs for $chroms[$chri] with:\n$fullCmd\n";
    }
    else {
	$pm->start && next;
	my $now = strftime("%F %T", localtime);
	warn "I: $now - $0 GATK4-GenotypeGVCFs for $chroms[$chri] starting\n";
	my $retVal = system($fullCmd);
	$now = strftime("%F %T", localtime);
	if ($retVal) {
	    warn "E: $now - $0 GATK4-GenotypeGVCFs for $chroms[$chri] FAILED ($?)!  INSPECT THE LOGFILE $log\n";
	}
	else {
	    warn "I: $now - $0 GATK4-GenotypeGVCFs for $chroms[$chri] completed successfully\n";
	}
        $pm->finish;
    }
}
$pm->wait_all_children;


#############################################
# merge chrom-specific VCFs with $bcftools, producing $outFile

# $outFile has been stripped of .vcf and .gz if present, add them back
$outFile .= ".vcf.gz";

$cmd = "$bcftools concat --output $outFile --output-type z";

# for speed: multiple compression threads
$cmd .= " --threads $jobs";

# if it's too slow I could also try --naive

# list of temp chrom-specific files, in correct order
$cmd .= " ".join(' ', @tempVcfs);

if (! $real){
    warn "I $0: dryrun, would run bcftools-concat with:\n$cmd\n";
}
else {
    my $now = strftime("%F %T", localtime);
    warn "I: $now - $0 bcftools-concat starting\n";
    my $retVal = system($cmd);
    $now = strftime("%F %T", localtime);
    if ($retVal) {
	# damn, now you wish we didn't use CLEANUP...
	die "E: $now - $0 bcftools-concat FAILED ($?)! Full command-line was:\n$cmd\n";
    }
    else {
	warn "I: $now - $0 bcftools-concat completed successfully\n";
	warn "I: $now - $0 ALL DONE\n";
    }
}
