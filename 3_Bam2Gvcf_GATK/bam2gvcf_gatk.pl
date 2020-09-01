#!/usr/bin/perl


# 17/05/2019
# NTM

# 13/06/2020
# OB

# Calls variants on BAM files and produces corresponding GVCF files using GATK4.
# The body of this script is quite similar to that of bam2gvcf_strelka.pl

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Parallel::ForkManager;



#############################################
## options / params from the command-line

my $USAGE = '
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--outdir string [required] : dir where GVCF files will be created
--indir string [defaults to a working path on luxor] : 
    subdir containing the BAMs
--first int [default = 50] : first grexomeNum to process (>= 50)
--last int [default = 9999] : last grexomeNum to process (>= $first)
--gatk [default path in luxor] : full path to gatk executable
--genome string [defaults to searching in NTM\'s standard places] : path+filename
    of the reference genome in fasta format, subdir must also contain chroms-only BED
--jobs N [default = 15] : number of samples to process at the same time. It is NOT 
    recommended to go below 4, because the processing time will decrease significantly
    (on luxor : 15 jobs -> ~20min/sample; 1 job -> ~2hours/sample). The parameter should 
    be adjusted depending on CPU capabilties, default value = 15 was tweaked on luxor.
--real : actually do the work, otherwise this is a dry run, just print 
    info on what would be done
--disable-log [default=false] : disable the creation of a log file (one per sample).
--help : print this USAGE';



# dir where GVCF-containing subdirs will be created, no default => MUST BE PROVIDED
my $outDir = "";

# subdir where BAMs are stored, checking for defaults in luxor and dahu
my $inDir = "";

# first and last grexomeNums to process, default to everything!
my ($firstGrex, $lastGrex) = (50,9999);

# full path for GATK executable, checking for defaults in luxor and dahu
my $gatk = "";

# full path to ref genome, must be the one used for producing the BAMs
my $refGenome;

# number of parallel jobs to run, get from command-line --jobs, default to 15
my $jobs = 15;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# disable log option, default to false
my $disable_log = '';

# help: if true just print $USAGE and exit
my $help = '';

GetOptions ("outdir=s" => \$outDir,
	    "indir=s" => \$inDir,
	    "first=i" => \$firstGrex,
	    "last=i" => \$lastGrex,
	    "gatk=s" => \$gatk,
	    "genome=s" => \$refGenome,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "disable-log" => \$disable_log,
	    "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) &&
    die "$USAGE\n\n";
($outDir) || 
    die "$USAGE\n\nE: you MUST specify the dir where GVCF-containing subdirs will be created, with --outdir\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E: outDir $outDir doesn't exist as a dir but can't be created\n";
if ($inDir){
    (-d $inDir || die "E : inDir specified is not a folder!");
}
else{
    if (-e "/data/nthierry/PierreRay/BAMs_All_Selected/") {
        $inDir = "/data/nthierry/PierreRay/BAMs_All_Selected/";
    }
    elsif (-e "/bettik/nthierry/BAMs_All_Selected/") {
        $inDir = "/bettik/nthierry/BAMs_All_Selected/";
    }
    else {
        die "E: can't find BAMs in default paths, use --indir option";
    }
} 
    
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


# make sure gatk executable exists, searching in default paths for luxor and dahu
if ($gatk) {
    (-f "$gatk") ||
	die "E: the provided gatk executable doesn't exist or can't be read\n";
}
elsif (-f "/usr/local/bin/gatk"){
#elsif (!(`which gatk`)){ # gatk is found in $PATH
    $gatk = "gatk";
}
elsif (-f "/home/bencheko/Software/latest-gatk/gatk") { # default path in dahu
    $gatk = "/home/bencheko/Software/latest-gatk/gatk";
}
else {
    die "E: cannot find gatk executable in the default paths, you must use eg --gatk path/to/ref/gatk";
}


#tmp dir
my $tmpDir = " ";
if (-d "/home/nthierry/PierreRay_DATA/RunSecondaryAnalyses/RamDisk"){
    $tmpDir = " --tmp-dir /home/nthierry/PierreRay_DATA/RunSecondaryAnalyses/RamDisk ";
}

#############################################
## generating GQ bands and corresponding string

my $gqb_option = "";
my $cpt = 5;
while ($cpt<30){
    $gqb_option .= "-GQB $cpt ";
    $cpt += 3;
}
while ($cpt<50){
    $gqb_option = $gqb_option."-GQB $cpt ";
    $cpt = $cpt + 5;
}
while ($cpt<80){
    $gqb_option = $gqb_option."-GQB $cpt ";
    $cpt = $cpt + 10;
}




#############################################
## The following log files are stored in $outDir/log :
#  done.log : in append mode, records every success ever
#  failed_$time : records failures, log file unique per run to remove useless logs quickly
(-e "$outDir/log") || mkdir("$outDir/log") || 
    die "E: folder $outDir/log can't be created\n";

my $done_log = "$outDir/log/done.log";
open(DONE, ">> $done_log") ||
    die("E : Can't open nor create $done_log");
my $now = strftime("%Y-%m-%d-%H-%M", localtime);
my $fail_log = "$outDir/log/fail_$now.log";
open(FAIL, ">> $fail_log") ||
    die("E : Can't open nor create $fail_log");



#############################################
## preparing the list of commands (one per sample)

$now = strftime("%F %T", localtime);
print "I: $now - $0 STARTING TO WORK\n";

my $manager = new Parallel::ForkManager($jobs);
foreach (my $gNum = $firstGrex; $gNum<=$lastGrex; $gNum++) {
    my $cmd = ""; # GATK commandline for a single grexome
	  
    # grexome name
	my $grexome = $gNum;

    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";

    # make sure we have bam and bai files for $grexome, otherwise skip
    my $bam = "$inDir/${grexome}.bam";
    ((-e $bam) && ((-e "$bam.bai") || (-e "$inDir/${grexome}.bai"))) || ((print "W: $grexome was asked for but we don't have a BAM and BAI for it in $inDir, skipping this grexome\n") && next);

    # gvcf to produce
    my $gvcf = "$outDir/${grexome}.g.vcf.gz";
    # avoid overwriting a good file!
    (-e "$outDir/${grexome}.g.vcf.gz") && (print "W: $grexome was asked for but there is already a ${grexome}.g.vcf.gz. Please delete it if needed and rerun for this exome.\n") && next;
    
    # TODO: compare runtimes with larger -Xmx ?
    # remove StrandBiasBySample annotation (useless and pollutes the logs)
    $cmd = "$gatk --java-options \"-Xmx4g\" HaplotypeCaller -R $refGenome -I $bam -O $gvcf --emit-ref-confidence GVCF --native-pair-hmm-threads 4 --verbosity INFO $tmpDir ";
    # maybe add "--verbosity WARNING" (at least in prod)

    my $sleep = ($gNum - 50)*0; ## experimental..

    # adding gqb option
    $cmd = $cmd.$gqb_option;
    #$cmd = "sleep $sleep && ".$cmd.$gqb_option; ## also experimental

    # log option
    if (! $disable_log){
        my $log = "$outDir/${grexome}.log";
        $cmd = $cmd." &>$log";
    }
    else{
        $cmd = $cmd." &> /dev/null";
    }

    # dryrun if not real
    if (! $real) {
        print "I: dryrun, would run gatk for $grexome with: $cmd\n";
    }
    else {
        $manager->start and next;
        if (system($cmd) != 0) {
            $now = strftime("%F %T", localtime);
            print FAIL "\n$now \n\n $cmd failed.\n";
        }else{
            $now = strftime("%F %T", localtime);
            print DONE "\n$now \n\n $cmd done.\n";
        }
        $manager->finish;
    }
}
$manager->wait_all_children;

close FAIL;
close DONE;

$now = strftime("%F %T", localtime);
print "I: $now - $0 ALL DONE\n";


