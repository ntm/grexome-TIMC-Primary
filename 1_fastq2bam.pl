#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2025
#
# This file is part of grexome-TIMC-Primary, written by Nicolas Thierry-Mieg
# (CNRS, France) Nicolas.Thierry-Mieg@univ-grenoble-alpes.fr
#
# This program is free software: you can redistribute it and/or modify it under 
# the terms of the GNU General Public License as published by the Free Software 
# Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY;
# without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program.
# If not, see <https://www.gnu.org/licenses/>.
############################################################################################


# 28/05/2019
# NTM

# Process a bunch of FASTQ paired-end files:
# - trim adaptors and filter low-quality reads (with fastp)
# - align reads (with bwa-mem2 if available, bwa otherwise)
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
# In particular, for GRCh38 we use $bwakitPostalt, which produces a bunch
#  of *hla* files, we don't use them but we keep them anyways: they're not 
# huge and and will be there if we ever want to do HLA typing.
# For this reason, when working on GRCh38 the reference genome should be
# produced by run-gen-ref from bwa-kit.
#
# Note on error-handling:
# - all logging goes to stderr;
# - we die for blatant problems in the prep stage;
# - if prep was OK and we started processing samples, we log E: messages but
#   never die. This avoids aborting the whole job (eg on a cluster) when
#   we have just a few samples failing. At the end we log the number of
#   errors and warnings that occurred.
#
# See $USAGE

use strict;
use warnings;
use File::Basename qw(basename);
use File::Path qw(remove_tree);
use Getopt::Long;
use POSIX qw(strftime);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## hard-coded stuff that shouldn't change

# programs used
my $fastp = "fastp";
my $samblaster = "samblaster";
my $samtools = "samtools";
# for BWA we'll use the first element of @bwas that exists (ie currently bwa-mem2
# if we find it, bwa otherwise)
my @bwas = ("bwa-mem2", "bwa");

# require bash so we can use -o pipefail
my $bash = "bash";

#############################################
## options / params from the command-line


# subdir where FASTQs can be found (required)
my $inDir = '';

# comma-separated list of samples (FASTQs) to process (required)
my $samples = '';

# subdir where BAMS will be created (required)
my $outDir = '';

# tmpDir, must not pre-exist and will be rm'd, faster is better
# (ramdisk or at least SSD)
my $tmpDir;

# path to programs used (fastp etc...), empty string seaches in $PATH,
# otherwise it must be a slash-terminated path (but we add trailing 
# slash if needed)
# If your binaries are in various dirs all over the place, you should
# create a subdir somewhere and symlink all the required binaries there,
# then use that subdir as binPath
my $binPath = '';

# also need path to bwa-kit subdir (with k8 and bwa-postalt.js),
# gets its own variable because it should never be in PATH.
# Leave empty to ignore the bwakit-Postalt step (eg non-human data)
my $bwakit = '';

# path+filename of ref genome, currently we recommend the full GRCh38 with
# decoy+alts+unmapped, as produced by Heng Li's run-gen-ref (from bwa-kit)
my $genome;

# max number of threads
my $numThreads = 4;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';


my $USAGE = "\nProcess a bunch of FASTQ paired-end files:
1. trim adaptors and filter low-quality reads (with fastp)
2. align reads (with bwa-mem2 if available, bwa otherwise)
3. mark dupes (with samblaster)
4. sort the BAM (with samtools)

Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir : subdir containing the FASTQs
--samples : comma-separated list of sampleIDs to process, for each sample there should be 
          a pair of FASTQ files in indir called [sample]_1.fq.gz and [sample]_2.fq.gz
--outdir : subdir where BAMs and accessory files will be created
--tmpdir : subdir where tmp files will be created, must not pre-exist and will be removed after execution
--binpath : path where binaries $fastp, $samblaster, $samtools and 
    ".join('/',@bwas)." can be found, leave empty to search in PATH
--bwakit : when aligning on GRCh38, path where k8 and bwa-postalt.js (from bwa-kit) can
          be found; if not provided, ignore bwakit-PostAlt (eg with non-human data)
--genome : ref genome fasta, with path, must be indexed with 'bwa-mem2 index' and/or 'bwa index'
--threads N [default = $numThreads] : number of threads for BWA, and also for fastp and samtools if <= 16
    (but if > 16 fastp and samtools uses only 16 threads)
--real : actually do the work, otherwise this is a dry run, just print info on what would be done
--help : print this USAGE";

GetOptions ("indir=s" => \$inDir,
            "samples=s" => \$samples,
            "outdir=s" => \$outDir,
            "tmpdir=s" => \$tmpDir,
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

($tmpDir) || die "E $0: you must provide a tmpDir\n";
(-e $tmpDir) &&
    die "E $0: found argument $tmpDir but it already exists, remove it or choose another name.\n";
mkdir($tmpDir) || die "E $0: cannot mkdir tmpDir $tmpDir\n";

# slash-terminate $binPath if it's not empty
($binPath) && (($binPath  =~ m~/$~)  || ($binPath .= "/"));

# make sure ref genome exists
($genome) || die "E $0: you must provide a ref genome fasta file\n";
(-f $genome) || die "E $0: provided genome fasta file doesn't exist\n";

# actual bwa-postalt command (use k8 to interpret the js)
my $bwakitPostalt;
if ($bwakit) {
    $bwakitPostalt = "$bwakit/k8 $bwakit/bwa-postalt.js";
    (`$bwakitPostalt -v` =~ /^r\d+$/) ||
        die "E $0: bwakitPostalt test doesn't run as expected, maybe fix bwakitPath in config.pm, command run: $bwakitPostalt -v\n";
    (-f "$genome.alt") ||
        die "E $0: provided ref genome found but we also need $genome.alt for bwa-postalt, ".
        "as produced by Heng Li's run-gen-ref (from bwa-kit)\n";
}

# make sure all progs can be found
system("which $bash &> /dev/null") && die "E $0: the bash executable $bash can't be found\n";
system("which $binPath$fastp &> /dev/null") && die "E $0: the fastp executable $fastp can't be found\n";
system("which $binPath$samblaster &> /dev/null") && die "E $0: the samblaster executable $samblaster can't be found\n";
system("which $binPath$samtools &> /dev/null") && die "E $0: the samtools executable $samtools can't be found\n";

my $bwa = "";
foreach my $b (@bwas) {
    if (system("which $binPath$b &> /dev/null") == 0) {
        # make sure genome is indexed for this flavor of BWA
        if ($b =~ /mem/) {
            # bwa-mem2
            if ((-f "$genome.bwt.2bit.64") && (-f "$genome.0123") && (-f "$genome.pac") && 
                (-f "$genome.ann") && (-f "$genome.amb")) {
                $bwa = $b;
                last;
            }
            else {
                warn "W $0: found $b but ref genome $genome isn't indexed for it, please use '$b index'\n";
                next;
            }
        }
        else {
            # classic BWA index
            if ((-f "$genome.bwt") && (-f "$genome.sa") && (-f "$genome.pac") && 
                (-f "$genome.ann") && (-f "$genome.amb")) {
                $bwa = $b;
                last;
            }
            else {
                warn "W $0: found $b but ref genome $genome isn't indexed for it, please use '$b index'\n";
                next;
            }
        }
    }
}
($bwa) || die "E $0: cannot find any usable (indexed genome) BWA executable among (".join(',',@bwas).")\n";

# ok, prepend binPath
$fastp = "$binPath$fastp";
$bwa = "$binPath$bwa";
$samblaster = "$binPath$samblaster";
$samtools = "$binPath$samtools";

# number of threads for samtools and fastp: capped at 16
my $numThreadsCapped = 16;
($numThreads < $numThreadsCapped) && ($numThreadsCapped = $numThreads);

# number of samples for which we got errors (resp warnings)
my $nbErrors = 0;
my $nbWarnings = 0;


#############################################
## build list of sanity-checked samples to process
# key == sampleID, value == 1
my %samples;
# single timestamp for this prep stage, should all happen in 1s
my $now = strftime("%F %T", localtime);
foreach my $sample (split(/,/, $samples)) {
    if ($samples{$sample}) {
        warn "W $now: $0 - sample $sample was specified twice, is that a typo? Ignoring the dupe\n";
        $nbWarnings++;
        next;
    }
    # fastq files, this MUST MATCH $f1 and $f2 that we process later
    my $f1 = "$inDir/${sample}_1.fq.gz";
    my $f2 = "$inDir/${sample}_2.fq.gz";
    if ((! -f $f1) || (! -f $f2)) {
        warn "W $now: $0 - sample $sample was specified but we don't have a pair of FASTQs for it in $inDir, skipping\n";
        $nbWarnings++;
        next;
    }
    my $bam = "$outDir/$sample.bam";
    if (-e $bam) {
        warn "W $now: $0 - sample $sample was specified but we already have the BAM $bam, remove it to re-process this sample, skipping\n";
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
    # discard json output, keep HTML output (detailed), and log stderr
    # other stuff is left at default, ie: no quality trimming, quality 
    # filtering filters reads with >5 N's or >40% low-qual (Q<15) bases,
    # length filtering filters reads shorter than 15 bp
    my $com = "$fastp --stdout --in1 $f1 --in2 $f2 --detect_adapter_for_pe --json /dev/null";
    $com .= " --html ${outFile}_fastp.html --thread $numThreadsCapped 2> ${outFile}_fastp.log | ";
    
    # BWA: -p (interleaved fastq), -R to add read group info,
    # -K 100000000 to make bwa reproducible (otherwise you can get different 
    # results when running with different numbers of threads!)
    $com .= "$bwa mem -p -t$numThreads -R \'\@RG\\tID:$sample\\tSM:$sample\' -K 100000000 $genome - 2> ${outFile}_bwa.log | ";

    # samblaster: nothing special
    $com .= "$samblaster 2> ${outFile}_samblaster.log |";

    # bwa-kit run-bwamem has a step for dealing correctly with ALT contigs (bwa-postalt.js),
    # we run that script too if --bwakit was provided
    # (see https://github.com/lh3/bwa/blob/master/README-alt.md )
    ($bwakitPostalt) && ($com .= "$bwakitPostalt -p $outFile.hla $genome.alt |");
    
    # sort with samtools
    $com .= "$samtools sort -\@ $numThreadsCapped -m1G -T $tmpDir -o $outFile.bam - ";

    # with bash -o pipefail, the return status of $com will be non-zero if any
    # piped component of $com fails
    # -> need to protect double-quotes
    $com =~ s/"/\\"/g;
    $com = "$bash -o pipefail -c \" $com \"";
    
    $now = strftime("%F %T", localtime);
    if (! $real) {
        warn "I $now: $0 - dryrun, would run: $com\n";
    }
    else {
        warn "I $now: $0 - starting processing of $sample with command: $com\n";
        if (system($com)) {
            # non-zero exit status
            $now = strftime("%F %T", localtime);
            warn "E $now: $0 - processing of $sample exited with non-zero status. Something went wrong, investigate!\n";
            $nbErrors++;
            # remove (corrupt) bamfile if it exists
            unlink("$outFile.bam");
            next;
        }

        $now = strftime("%F %T", localtime);
        warn "I $now: $0 - done aligning $sample, indexing\n";
        if(system("$samtools index $outFile.bam")) {
            # non-zero exit status
            $now = strftime("%F %T", localtime);
            warn "E $now: $0 - samtools index $sample exited with non-zero status. Strange because bam production succeeded, investigate!\n";
            $nbErrors++;
        }
    }
}

remove_tree($tmpDir);
(-e $tmpDir) && 
    warn "E $now: $0 - all done but cannot rmdir tmpDir $tmpDir, why?\n";

$now = strftime("%F %T", localtime);
if ($nbErrors) {
    warn "E $now: $0 - finished but $nbErrors ERRORS DETECTED, I was running ".join(" ", $0, @ARGV)."\n";
    exit(1);
}
elsif ($nbWarnings) {
    warn "W $now: $0 - finished but $nbWarnings WARNINGS need verification, I was running ".join(" ", $0, @ARGV)."\n";
}
else {
    warn "I $now: $0 - finished SUCCESSFULLY, I was running ".join(" ", $0, @ARGV)."\n";
}
