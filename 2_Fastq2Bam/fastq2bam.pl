#!/usr/bin/perl


# 28/05/2019
# NTM

# Process a bunch of grexome FASTQ paired-end files:
# - trim adaptors and filter low-quality reads (with fastp)
# - align reads (with bwa-mem)
# - mark dupes (with samblaster)
# - sort the BAM (with samtools)
#
# We are quite stringent on the naming convention for the FASTQ files:
# they must be PE and MUST BE CALLED grexome\d\d\d\d_[12].fq.gz .
# BAM files will be created in $outDir (skipping any sample with
# a pre-existing BAM).
#
# This is inspired by bwa-kit.
# In particular we use $bwakitPostalt which produces a bunch of *hla* 
# files, we don't use them but we keep them anyways: they're not 
# huge and and will be there if we ever want to do HLA typing
#
# We log to stderr and die for blatant problems in the prep stage;
# if prep was OK and we started processing samples, we log to stdout 
# and never die, check stdout for E: and W: messages.
#
# See $USAGE

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);


#############################################
## hard-coded stuff that shouldn't change

# programs used
my $fastp = "fastp";
my $bwa = "bwa";
my $samblaster = "samblaster";
my $samtools = "samtools";


#############################################
## options / params from the command-line

my $USAGE = '
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--outdir string [required] : subdir where BAMs will be created
--indir string [defaults to a working path on luxor and fauve] : subdir 
    containing the grexomized FASTQs
--first int [default = 50] : first grexomeNum to process (>= 50)
--last int [default = 9999] : last grexomeNum to process (>= $first)
--binpath string [default ""]: path where binaries fastp etc can be found,
    leave empty to search in PATH
--bwakit string [default "~/Software/BWA-kit/bwa.kit/"] : path where k8 and
    bwa-postalt.js (from bwa-kit) can be found
--genome string [defaults to searching in NTM\'s standard places] : path+filename
    of the reference genome in fasta format, must be indexed with "bwa index"
--threads N [default = 4] : number of threads for BWA, and also for fastp
    and samtools if <= 4 (but if > 4 fastp and samtools use only 4 threads)
--real : actually do the work, otherwise this is a dry run, just print 
    info on what would be done
--help : print this USAGE';


# subdir where BAMS will be created, no default => MUST BE PROVIDED as arg
my $outDir = "";

# subdir where grexomized FASTQs are stored, defaults to current
# correct dir on luxor and also on fauve (and others)
my $inDir = "/data/nthierry/PierreRay/FASTQs_All_Grexomized/";
(-e "/home/nthierry/sshMounts/luxor/") &&
    ($inDir = "/home/nthierry/sshMounts/luxor/$inDir");

# first and last grexomeNums to process, default to everything!
my ($firstGrex, $lastGrex) = (50,9999);

# path to programs used (fastp etc...), empty string seaches in $PATH,
# otherwise it must be a slash-terminated path (but we add trailing 
# slash if needed)
my $binPath = "";

# also need path to bwa-kit subdir (with k8 and bwa-postalt.js),
# gets its own variable because it should never be in PATH
# set a default that works on fauve, luxor, krakenator
my $bwakit = "~/Software/BWA-kit/bwa.kit/";

# path+filename of ref genome, currently this should be the full GRCh38 with
# decoy+alts+unmapped, as produced by Heng Li's run-gen-ref (from bwa-kit)
my $refGenome;

# number of threads, get from command-line -t, default to 4
my $numThreads = 4;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

GetOptions ("outdir=s" => \$outDir,
	    "indir=s" => \$inDir,
	    "first=i" => \$firstGrex,
	    "last=i" => \$lastGrex,
	    "binpath=s" => \$binPath,
	    "bwakit=s" => \$bwakit,
	    "genome=s" => \$refGenome,
	    "threads=i" => \$numThreads, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) &&
    die "$USAGE\n\n";
($outDir) || 
    die "$USAGE\n\nE: you MUST specify the dir where BAMs will be created, with --outdir\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E: outDir $outDir doesn't exist as a dir but can't be created\n";
(-d $inDir) || 
    die "E: inDir $inDir doesn't exist\n";
(($firstGrex >= 50) && ($lastGrex >= $firstGrex)) ||
    die "E: first grexomeNum must be >=50 and last must be >=first\n";

# bring $lastGrex down to the largest existing grexome*_1.fq.gz
# (mostly in case we have the default 9999)
while($lastGrex > $firstGrex) {
    my $grexome = $lastGrex;
    # left-pad with zeroes to 4 digits
    ($grexome < 10) && ($grexome = "0$grexome");
    ($grexome < 100) && ($grexome = "0$grexome");
    ($grexome < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";
    (-f "$inDir/${grexome}_1.fq.gz") && last;
    $lastGrex--;
}

# slash-terminate $binPath if it's not empty
($binPath) && (($binPath  =~ m~/$~)  || ($binPath .= "/"));

# make sure all progs can be found
(`which $binPath$fastp` =~ /$fastp$/) || die "E: the fastp executable $fastp can't be found\n";
(`which $binPath$bwa` =~ /$bwa$/) || die "E: the bwa executable $bwa can't be found\n";
(`which $binPath$samblaster` =~ /$samblaster$/) || die "E: the samblaster executable $samblaster can't be found\n";
(`which $binPath$samtools` =~ /$samtools$/) || die "E: the samtools executable $samtools can't be found\n";
# ok, prepend binPath
$fastp = "$binPath$fastp";
$bwa = "$binPath$bwa";
$samblaster = "$binPath$samblaster";
$samtools = "$binPath$samtools";

# actual bwa-postalt command (use k8 to interpret the js)
my $bwakitPostalt = "$bwakit/k8 $bwakit/bwa-postalt.js";
(`$bwakitPostalt -v` =~ /^r\d+$/) ||
    die "E: bwakitPostalt test doesn't run as expected, maybe fix bwakit subdir, command run: $bwakitPostalt -v\n";

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

# here we know $refGenome exists, but BWA needs it indexed
((-f "$refGenome.bwt") && (-f "$refGenome.pac") && (-f "$refGenome.sa") && 
 (-f "$refGenome.ann") && (-f "$refGenome.amb")) ||
    die "E: the reference genome $refGenome isn't indexed, please use 'bwa index'\n";


# number of threads for fastp and samtools: capped at 4
my $numThreadsCapped = 4;
($numThreads < $numThreadsCapped) && ($numThreadsCapped = $numThreads);


#############################################
## process each sample of interest

# number of samples for which we got errors (resp warnings)
my $nbErrors = 0;
my $nbWarnings = 0;

foreach my $gNum ($firstGrex..$lastGrex) {
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";

    # sanity: make sure we have fastq files for $grexome
    my $f1 = "$inDir/${grexome}_1.fq.gz";
    my $f2 = "$inDir/${grexome}_2.fq.gz";
    if ((! -f $f1) || (! -f $f2)) {
	print "W: $grexome was asked for but we don't have a pair of FASTQs for it in $inDir, skipping this grexome\n";
	$nbWarnings++;
	next;
    }

    # will create $outFile* files, mainly .bam but also 
    # some log files starting with $outFile
    my $outFile = "$outDir/$grexome";
    
    # if $outFile* exists we skip this grexome with a warning
    my @files = glob($outFile."*");
    if (@files > 0) {
	print "W: outfiles exist for grexome $grexome in outdir $outDir, remove them to process this sample, skipping\n";
	$nbWarnings++;
	next;
    }

    # fastp: enable autodetection of adaptors (in addition to overlap analysis),
    # discard json output, keep HTML output (detailed) and log stderr
    # other stuff is left at default, ie: no quality trimming, quality 
    # filtering filters reads with >5 N's or >40% low-qual (Q<15) bases,
    # length filtering filters reads shorter than 15 bp
    my $com = "$fastp --stdout --in1 $f1 --in2 $f2 --detect_adapter_for_pe --json /dev/null --html ${outFile}_fastp.html --thread $numThreadsCapped 2> ${outFile}_fastp.log | ";
    
    # BWA: -p (interleaved fastq), -R to add read group info,
    # -K 100000000 to make bwa reproducible (otherwise you can get different 
    # results when running with different numbers of threads!)
    $com .= "$bwa mem -p -t$numThreads -R \'\@RG\\tID:$grexome\\tSM:$grexome\' -K 100000000 $refGenome - 2> ${outFile}_bwa.log | ";

    # samblaster: nothing special
    $com .= "$samblaster 2> ${outFile}_samblaster.log |";

    # bwa-kit run-bwamem has a step for dealing correctly with ALT
    # contigs (bwa-postalt.js), we run that script too
    # (see https://github.com/lh3/bwa/blob/master/README-alt.md )
    $com .= "$bwakitPostalt -p $outFile.hla $refGenome.alt |";
    
    # sort with samtools
    $com .= "$samtools sort -\@ $numThreadsCapped -m1G -o $outFile.bam - ";

    if (! $real) {
	print "I: dryrun, would run: $com\n";
    }
    else {
	my $now = strftime("%F %T", localtime);
	print "I: $now - starting processing of $grexome with command: $com\n";
	if (system($com)) {
	    # non-zero exit status
	    print "E: processing of $grexome exited with non-zero status. Something went wrong, investigate!\n";
	    $nbErrors++;
	    # don't even try to index the bam
	    next;
	}

	$now = strftime("%F %T", localtime);
	print "I: $now - done aligning $grexome, indexing\n";
	if(system("$samtools index $outFile.bam")) {
	    # non-zero exit status
	    print "E: samtools index $grexome exited with non-zero status. Strange because processing didn't seem to fail, investigate!\n";
	    $nbErrors++;
	}
    }
}

my $now = strftime("%F %T", localtime);
if ($nbErrors) {
    print "E: $now - finished but $nbErrors ERRORS DETECTED, I was running ".join(" ", $0, @ARGV)."\n";
}
elsif ($nbWarnings) {
    print "W: $now - finished but $nbWarnings WARNINGS need verification, I was running ".join(" ", $0, @ARGV)."\n";
}
else {
    print "I: $now - finished SUCCESSFULLY, I was running ".join(" ", $0, @ARGV)."\n";
}
