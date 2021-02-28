#!/usr/bin/perl


# 28/05/2019
# NTM

# Process a bunch of FASTQ paired-end files:
# - trim adaptors and filter low-quality reads (with fastp)
# - align reads (with bwa-mem)
# - mark dupes (with samblaster)
# - sort the BAM (with samtools)
#
# We are quite stringent on the naming convention for the FASTQ files:
# for each listed $sample, we expect a single pair of FASTQ files living
# in $indir called ${sample}_1.fq.gz and ${sample}_2.fq.gz .
# BAM files will be created in $outDir (skipping any sample with
# a pre-existing BAM).
#
# This is inspired by bwa-kit.
# In particular we use $bwakitPostalt, which produces a bunch of *hla* 
# files, we don't use them but we keep them anyways: they're not 
# huge and and will be there if we ever want to do HLA typing.
# For this reason the reference genome should be produced by
# run-gen-ref from bwa-kit.
#
# We log to stderr and die for blatant problems in the prep stage;
# if prep was OK and we started processing samples, we log to stdout 
# and never die, check stdout for E: and W: messages.
# This was done to avoid aborting the whole job (eg on a cluster) when
# we have just a few samples failing.
#
# See $USAGE

use strict;
use warnings;
use File::Basename qw(basename);
use Getopt::Long;
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## hard-coded stuff that shouldn't change

# programs used
my $fastp = "fastp";
my $bwa = "bwa";
my $samblaster = "samblaster";
my $samtools = "samtools";


#############################################
## options / params from the command-line


# subdir where FASTQs can be found (required)
my $inDir = '';

# comma-separated list of samples (FASTQs) to process (required)
my $samples = '';

# subdir where BAMS will be created (required)
my $outDir = '';

# path to programs used (fastp etc...), empty string seaches in $PATH,
# otherwise it must be a slash-terminated path (but we add trailing 
# slash if needed)
my $binPath = '';

# also need path to bwa-kit subdir (with k8 and bwa-postalt.js),
# gets its own variable because it should never be in PATH
# set a default that works on fauve, luxor, krakenator
my $bwakit = "~/Software/BWA-kit/bwa.kit/";

# path+filename of ref genome, currently we recommend the full GRCh38 with
# decoy+alts+unmapped, as produced by Heng Li's run-gen-ref (from bwa-kit)
my $genome;

# number of threads, get from command-line -t, default to 4
my $numThreads = 4;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nProcess a bunch of FASTQ paired-end files:
1. trim adaptors and filter low-quality reads (with fastp)
2. align reads (with bwa-mem)
3. mark dupes (with samblaster)
4. sort the BAM (with samtools)

Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir string : subdir containing the FASTQs
--samples : comma-separated list of sampleIDs to process, for each sample there should be 
	  a pair of FASTQ files in indir called [sample]_1.fq.gz and [sample]_2.fq.gz
--outdir string : subdir where BAMs and accessory files will be created
--binpath string [default '']: path where binaries fastp etc can be found,
    leave empty to search in PATH
--bwakit string [default '~/Software/BWA-kit/bwa.kit/'] : path where k8 and
    bwa-postalt.js (from bwa-kit) can be found
--genome string [no default] : ref genome fasta, with path, must be indexed with 'bwa index'
--threads N [default = 4] : number of threads for BWA, and also for fastp and
    samtools if <= 4 (but if > 4 fastp and samtools use only 4 threads)
--real : actually do the work, otherwise this is a dry run, just print 
    info on what would be done
--help : print this USAGE";

GetOptions ("indir=s" => \$inDir,
	    "samples=s" => \$samples,
	    "outdir=s" => \$outDir,
	    "binpath=s" => \$binPath,
	    "bwakit=s" => \$bwakit,
	    "genome=s" => \$genome,
	    "threads=i" => \$numThreads, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($inDir) || 
    die "E $0: you MUST specify the dir where FASTQs can be found, with --indir\n$USAGE\n";
(-d $inDir) || 
    die "E $0: inDir $inDir doesn't exist\n";

($outDir) || 
    die "E $0: you MUST specify the dir where BAMs will be created, with --outdir\n$USAGE\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E $0: outDir $outDir doesn't exist as a dir and can't be created\n";

# slash-terminate $binPath if it's not empty
($binPath) && (($binPath  =~ m~/$~)  || ($binPath .= "/"));

# make sure all progs can be found
(`which $binPath$fastp` =~ /$fastp$/) || die "E $0: the fastp executable $fastp can't be found\n";
(`which $binPath$bwa` =~ /$bwa$/) || die "E $0: the bwa executable $bwa can't be found\n";
(`which $binPath$samblaster` =~ /$samblaster$/) || die "E $0: the samblaster executable $samblaster can't be found\n";
(`which $binPath$samtools` =~ /$samtools$/) || die "E $0: the samtools executable $samtools can't be found\n";
# ok, prepend binPath
$fastp = "$binPath$fastp";
$bwa = "$binPath$bwa";
$samblaster = "$binPath$samblaster";
$samtools = "$binPath$samtools";

# actual bwa-postalt command (use k8 to interpret the js)
my $bwakitPostalt = "$bwakit/k8 $bwakit/bwa-postalt.js";
(`$bwakitPostalt -v` =~ /^r\d+$/) ||
    die "E $0: bwakitPostalt test doesn't run as expected, maybe fix bwakit subdir, command run: $bwakitPostalt -v\n";

# make sure ref genome exists
($genome) || die "E $0: you must provide a ref genome fasta file\n";
(-f $genome) || die "E $0: provided genome fasta file doesn't exist\n";
(-f "$genome.alt") ||
    die "E $0: provided ref genome found but we also need $genome.alt for bwa-postalt, ".
    "as produced by Heng Li's run-gen-ref (from bwa-kit)\n";
# BWA needs the genome indexed
((-f "$genome.bwt") && (-f "$genome.pac") && (-f "$genome.sa") && 
 (-f "$genome.ann") && (-f "$genome.amb")) ||
    die "E $0: the reference genome $genome isn't indexed, please use 'bwa index'\n";
    
# number of threads for fastp and samtools: capped at 4
my $numThreadsCapped = 4;
($numThreads < $numThreadsCapped) && ($numThreadsCapped = $numThreads);

# number of samples for which we got errors (resp warnings)
my $nbErrors = 0;
my $nbWarnings = 0;


#############################################
## build list of sanity-checked samples to process
# key == sampleID, value == 1
my %samples;
foreach my $sample (split(/,/, $samples)) {
    if ($samples{$sample}) {
	print "W $0: sample $sample was specified twice, is that a typo? Ignoring the dupe\n";
	$nbWarnings++;
	next;
    }
    # fastq files, this MUST MATCH $f1 and $f2 that we process later
    my $f1 = "$inDir/${sample}_1.fq.gz";
    my $f2 = "$inDir/${sample}_2.fq.gz";
    if ((! -f $f1) || (! -f $f2)) {
	print "W $0: sample $sample was specified but we don't have a pair of FASTQs for it in $inDir, skipping\n";
	$nbWarnings++;
	next;
    }
    my $bam = "$outDir/$sample.bam";
    if (-e $bam) {
	print "W $0: sample $sample was specified but we already have the BAM $bam, remove it to re-process this sample, skipping\n";
	$nbWarnings++;
	next;
    }

    # AOK, sample will be processed
    $samples{$sample} = 1;
}


#############################################
## process each sample

foreach my $sample (sort keys(%samples)) {
    # fastq files, this MUST MATCH $f1 and $f2 sanity-checked above
    my $f1 = "$inDir/${sample}_1.fq.gz";
    my $f2 = "$inDir/${sample}_2.fq.gz";
    # we will create $outFile* files, mainly .bam but also 
    # some log files prefixed with $outFile
    my $outFile = "$outDir/$sample";

    # fastp: enable autodetection of adaptors (in addition to overlap analysis),
    # discard json output, keep HTML output (detailed) and log stderr
    # other stuff is left at default, ie: no quality trimming, quality 
    # filtering filters reads with >5 N's or >40% low-qual (Q<15) bases,
    # length filtering filters reads shorter than 15 bp
    my $com = "$fastp --stdout --in1 $f1 --in2 $f2 --detect_adapter_for_pe --json /dev/null --html ${outFile}_fastp.html --thread $numThreadsCapped 2> ${outFile}_fastp.log | ";
    
    # BWA: -p (interleaved fastq), -R to add read group info,
    # -K 100000000 to make bwa reproducible (otherwise you can get different 
    # results when running with different numbers of threads!)
    $com .= "$bwa mem -p -t$numThreads -R \'\@RG\\tID:$sample\\tSM:$sample\' -K 100000000 $genome - 2> ${outFile}_bwa.log | ";

    # samblaster: nothing special
    $com .= "$samblaster 2> ${outFile}_samblaster.log |";

    # bwa-kit run-bwamem has a step for dealing correctly with ALT
    # contigs (bwa-postalt.js), we run that script too
    # (see https://github.com/lh3/bwa/blob/master/README-alt.md )
    $com .= "$bwakitPostalt -p $outFile.hla $genome.alt |";
    
    # sort with samtools
    $com .= "$samtools sort -\@ $numThreadsCapped -m1G -o $outFile.bam - ";

    if (! $real) {
	print "I $0: dryrun, would run: $com\n";
    }
    else {
	my $now = strftime("%F %T", localtime);
	print "I $0: $now - starting processing of $sample with command: $com\n";
	if (system($com)) {
	    # non-zero exit status
	    print "E $0: processing of $sample exited with non-zero status. Something went wrong, investigate!\n";
	    $nbErrors++;
	    # don't even try to index the bam
	    next;
	}

	$now = strftime("%F %T", localtime);
	print "I $0: $now - done aligning $sample, indexing\n";
	if(system("$samtools index $outFile.bam")) {
	    # non-zero exit status
	    print "E $0: samtools index $sample exited with non-zero status. Strange because processing didn't seem to fail, investigate!\n";
	    $nbErrors++;
	}
    }
}

my $now = strftime("%F %T", localtime);
if ($nbErrors) {
    print "E $0: $now - finished but $nbErrors ERRORS DETECTED, I was running ".join(" ", $0, @ARGV)."\n";
}
elsif ($nbWarnings) {
    print "W $0: $now - finished but $nbWarnings WARNINGS need verification, I was running ".join(" ", $0, @ARGV)."\n";
}
else {
    print "I $0: $now - finished SUCCESSFULLY, I was running ".join(" ", $0, @ARGV)."\n";
}
