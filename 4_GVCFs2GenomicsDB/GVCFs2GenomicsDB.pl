#!/usr/bin/perl

# 16/10/2020
# NTM


# Import single-sample GATK GVCFs into a GenomicsDB: this is the
# intel-developed replacement for GATK's own CombineGVCFs (which
# never worked in our hands), and is therefore an alternative
# to our own mergeGVCFs.pl .
# This script creates a fresh DB if indir doesn't pre-exist,
# and updates an existing DB otherwise.


use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Basename qw(basename);
use FindBin qw($RealBin);
use File::Temp qw(tempdir);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded stuff that shouldn't change much

# RAM we give to GATK
my $ram = "192g";
# number of threads to use - RAM usage seems to scales linearly, and
# number of open file descriptors definitely scales linearly
my $threads = 12;

# NOTES:
# 1. GATK opens more than ($threads * $samples) files simultaneously,
#    when importing many samples you will probably need to increase
#    the "open files" limit in your shell to its hard limit (ulimit),
#    and then maybe decrease $threads until you're within the limit.
# 2. An alterative could be --batch-size (and maybe --consolidate).
# 3. I successfully imported 490 samples without batching using
#    $mem="288g" and $threads=8.



#############################################
## options / params from the command-line


# subdir where GVCFs can be found
my $inDir;

# subdir where GenomicsDB database will be created / updated
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

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = 'Import GVCFs into a GenomicsDB.
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir string [no default] : subdir containing the GVCFs
--outdir string [no default] : dir where GenomicsDB will be created/updated
--first int [$firstGrex] : first grexomeNum to process (>= 50)
--last int [$lastGrex] : last grexomeNum to process (>= $first)
--gatk [default to "gatk" which should be in PATH] : full path to gatk executable
--config string [$config] : your customized copy (with path) of the distributed *config.pm
--real : actually do the work, otherwise this is a dry run, just print info on what would be done
--help : print this USAGE';


GetOptions ("indir=s" => \$inDir,
	    "outdir=s" => \$outDir,
	    "first=i" => \$firstGrex,
	    "last=i" => \$lastGrex,
	    "gatk=s" => \$gatk,
	    "config=s" => \$config,
	    "real" => \$real,
	    "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenomeChromsBed fastTmpPath) );

($inDir) ||
    die "E $0: you MUST provide a dir where GVCFs can be found, with --indir\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

($outDir) || 
    die "E $0: you MUST specify the dir where the GenomicsDB will be created or updated, with --outdir\n";

(($firstGrex >= 50) && ($lastGrex >= $firstGrex)) ||
    die "E $0: first grexomeNum must be >=50 and last must be >=first\n";

# bring $lastGrex down to the largest existing grexome*.g.vcf.gz in inDir
# (mostly in case we have the default 9999)
while($lastGrex > $firstGrex) {
    my $grexome = $lastGrex;
    # left-pad with zeroes to 4 digits
    ($grexome < 10) && ($grexome = "0$grexome");
    ($grexome < 100) && ($grexome = "0$grexome");
    ($grexome < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";
    (-f "$inDir/${grexome}.g.vcf.gz") && last;
    $lastGrex--;
}

# make sure gatk executable is found, this test is disabled if
# we will be running GATK from a singularity container
($gatk =~ /singularity/) ||
    (`which $gatk` =~ /$gatk$/) ||
    die "E $0: cannot find 'gatk' (from GATK4 package), you must provide it with --gatk, you provided:\n$gatk\n";

# BED with chromosomes 1-22, X, Y, M
my $chromsBed = &refGenomeChromsBed();

# Create subdir in &fastTmpPath so we can CLEANUP when we
# are done, because GATK leaves a lot of sh*t behind
my $tmpDir = tempdir(DIR => &fastTmpPath(), CLEANUP => 1);

# test $outDir now to provide clear error messages and avoid calling GATK uselessly:
# if user wants to create a fresh DB $outDir cannot pre-exist,
# otherwise if $outDir exists it MUST contain a pre-existing genomicsDB
if (-d $outDir) {
    # callset.json should be present in any pre-existing GenomicsDB
    (-f "$outDir/callset.json") ||
	die "E $0: you provided a pre-existing outdir but it doesn't seem to be a GenomicsDB subdir.".
	" If you want to create a fresh DB outdir must not pre-exist\n";
}

my $now = strftime("%F %T", localtime);
my $logMess = "I: $now - $0 STARTING TO WORK, ";
if (-d $outDir) {
    $logMess .= "UPDATING PRE-EXISTING DB IN $outDir\n";
}
else {
    $logMess .= "CREATING NEW DB IN $outDir\n";
}
warn $logMess;

#############################################

# build GATK command-line
my $cmd = "$gatk  --java-options \"-Xmx$ram\" GenomicsDBImport";

if (-d $outDir) {
    $cmd .= " --genomicsdb-update-workspace-path $outDir";
}
else {
    # creating a fresh DB requires --intervals
    $cmd .= " --genomicsdb-workspace-path $outDir --intervals $chromsBed";
}

# use fast and big tmp storage
$cmd .= " --tmp-dir $tmpDir";

# max number of threads for import
$cmd .= " --max-num-intervals-to-import-in-parallel $threads";

# --batch-size 50  will help if we run out of memory or open file descriptors
# --consolidate is recommended if batch-size >= 100
# --reference ? why would the ref genome be needed?

foreach (my $gNum = $firstGrex; $gNum<=$lastGrex; $gNum++) {
    # grexome name
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";

    # make sure we have GVCF for $grexome, otherwise skip
    my $gvcf = "$inDir/${grexome}.g.vcf.gz";
    (-e $gvcf) ||
	((warn "W: no GVCF for $grexome in $inDir, skipping it\n") && next);
    
    $cmd .= " -V $gvcf";
}

if (! $real){
    warn "I: dryrun, would run GATK4-GenomicsDBImport with:\n$cmd\n";
}
else {
    my $retVal = system($cmd);
    $now = strftime("%F %T", localtime);
    if ($retVal) {
	die "E: $now - $0 running GATK4-GenomicsDBImport FAILED ($?)! Full command-line was:\n$cmd\n";
    }
    else {
	warn "I: $now - $0 ALL DONE\n";
    }
}

