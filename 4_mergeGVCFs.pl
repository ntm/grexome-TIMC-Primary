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


# 25/06/2019
# NTM

# Take as arg a filename containing a list of GVCF filenames (with full path, 
# possible bgzipped, one file per line) with one or more data columns (ie samples).
# These GVCFs should be produced by Strelka, GATK4, deepVariant or ElPrep, preferably
# cleaned up by filterBadCalls.pl, and/or by this program.
# If you feed random GVCF files you will probably need to adapt the code (although
# I try to be defensive): many aspects of the VCF spec are subject to interpretation...
# Produce to stdout a GVCF file, where:
# - header is copied from first file, except: 
#   all INFO descriptions except END and BLOCKAVG_min30p3a are stripped;
#   if --cleanheaders, ##contig headers are stripped except chroms \d+ and X,Y,M/MT (possibly
#      chr-prefixed);
#   ##mergeGVCFs= lines from all infiles are copied;
#   #CHROM line is modified by appending the identifiers of all samples.
# - any line whose FILTER column contains a key from %filtersApplied gets skipped.
# - lines whose FORMAT is "GT:PGT:PID:PS" are skipped (GATK4 produces these useless lines).
# - every variant gets normalized (but not left-aligned): 
#   * remove bases present at the end of REF and all ALTs;
#   * remove bases present at the start of REF and all ALTs, and increase POS accordingly;
#   * make sure there isn't another genotype call (other than non-var) at the new POS, 
#     and remove any non-var call that was there;
# - any remaining (adjusted) chrom:pos line from any infile will appear in the output, with:
#   ID set to '.';
#   the longest REF from all infiles;
#   all ALTs from all infiles, adjusted to fit the longest REF (append extra 
#      bases if needed);
#   QUAL and FILTER replaced by '.'
#   INFO replaced by '.' except for non-variant blocks:
#      * if all infiles have overlapping blocks, we produce a new block covering the intersection
#        (POS and END get updated, any other info such as BLOCKAVG_* from Strelka stays);
#      * otherwise we print individual lines, not a block;
#   FORMAT contains the union of FORMAT keys from all files (see below);
#   DATA columns contain the CORRECTED sample data for every sample.
#   CORRECTED means:
#     if a key was missing for a sample it gets '.'
#     GT gets corrected values for ALTs
#     GQ, GQX, DPI, DP, DPF, AF, SB, FT, PS, PGT, PID don't change
#     AD/ADF/ADR and VAF get 0 for new ALTs
#     PL gets correct new values, using 255 for missing alleles (see explanation in the code)
#
###########
# NOTES on FORMAT:
# From Strelka2 v2.9.10, for non-variant blocks it is always
# GT:GQX:DP:DPF:MIN_DP
# we also see this for non-variant sites (ie without END)
# the only formats for variant sites are:
# GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL
# GT:GQ:GQX:DPI:AD:ADF:ADR:FT:PL:PS
# GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL
# GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL:PS
#
# From GATK 4.1.8.1, in VCFs from GenotypeGVCFs via GenomicsDB we have only:
# GT:AD:DP:GQ:PGT:PID:PL:PS
# GT:AD:DP:GQ:PL
# in GATK GVCFs straight from HaplotypeCaller we have:
# GT
# GT:GQ:PL
# GT:PGT:PID:PS
# GT:AD:DP:GQ:PGT:PID:PL:PS:SB
# GT:AD:DP:GQ:PL:SB
# GT:DP:GQ:MIN_DP:PL
#
# From deepVariant 1.4.0 GVCFs we have only:
# GT:GQ:MIN_DP:PL (non-variant blocks)
# GT:GQ:DP:AD:VAF:PL (variant positions)
#
# AF is also accepted as a key (so input can be produced by filterBadCalls.pl).
#
# In the output files we can change the order of fields and
# create different combinations. But this script can read 
# it's own output.
###########
# NOTES on FILTERS
## -> Strelka2 v2.9.10:
## - frequent and seem all relevant: LowGQX LowDepth HighDPFRatio
##   UPDATE 30/10/2020 LowGQX is discarding variants that look good!
## - never seen without LowGQX or LowDepth: SiteConflict
## - rare, never seen with END, not the most relevant IMO: HighSNVSB
## - never seen in my data: IndelConflict NotGenotyped PloidyConflict
## -> GATK 4.1.8.1:
## - only available FILTER is LowQual, never seen in my data but should
##   be relevant
## -> deepVariant 1.4.0:
## - existing filters are RefCall (call should be homoref), LowQual (GQ should be low),
##   and NoCall (should be ./.), we'll just use LowQual since we already use it
###########
# NOTES on normalization and merging of variants:
# - as stated, I do NOT left-align variants, it's too much overhead. 
#   I am hoping that Strelka/GATK/DV do it correctly, or at least I'm hoping
#   they're consistent in their behavior. If they are consistent, merging
#   will work fine and variants can be left-aligned later if needed.
# - I also do NOT merge lines with different POS values, even if
#   the REFs overlap. I'm afraid doing this would result in large
#   REFs and huge lists of complex ALTs (most of which just describe
#   a SNV...).
# - if the input GVCFs went through filterBadCalls.pl the variants will
#   already be normalized, using the same code... However I don't want to
#   remove the normalization code from here, because I want mergeGVCFs.pl
#   to remain functional with vanilla GVCFs. Nevertheless, mergeGVCFs will
#   produce more useful results if the input GVCFs went through filterBadCalls
#   beforehand, since filterBadCalls does a lot of clean-up and compensates
#   for many variant-caller bugs^Hfeatures.


use strict;
use warnings;
use Getopt::Long;
use Parallel::ForkManager;
use POSIX qw(strftime);
use File::Basename qw(basename);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## hard-coded stuff that shouldn't change much

# filters to apply: any line whose FILTER value contains a key of %filtersApplied
# will be skipped.
# For performance reasons we do NOT check the FT fields for these keys, we assume
# that the infiles are single-sample GVCFs or were produced by this script
# with the same %filtersApplied.
# [see NOTES on FILTERS above for how I chose these]
my %filtersApplied = ('LowDepth'=>1, 'HighDPFRatio'=>1, 'LowQual'=>1);


# adaptive strategy parameters for $batchSize: without --batchSize, $batchSize will
# be adjusted in order to process ($jobs-1) batches in [$batchTimeLow,$batchTimeHigh] seconds
my $batchSizeAdaptive = 1;
# aim for between 5m and 30m for ($jobs-1) batches
my ($batchTimeLow, $batchTimeHigh) = (300,1800);


#############################################
## options / params from the command-line

# file containing list of GVCFs, no default
my $fileList;

# tmpDir, must not pre-exist and will be rm'd
my $tmpDir;

# number of parallel jobs to run
my $jobs = 16;

# $batchSize: max number of lines to read in a single batch from the first infile
# (number of lines from other infiles should be similar if the files have the same
# numbers of samples).
# Each batch is then processed by a worker thread.
# Increasing $batchSize increases the RAM consumption linearly, but should
# improve performance due to better read performance (because better OS buffering),
# less overhead from sub calls, process creations, tmpFiles creations, etc...
# We have used as low as 5k (merging 550 individual filtered exome-derived GVCFs),
# and as high as 100k (merging one multi-sample GVCF containing 540 samples with
# 10 individual GVCFs, all exome-derived).
# Providing it with --batchSize hard-sets $batchSize, could be useful if you lack RAM.
# Without --batchSize we start with a low value and use an adaptive strategy
my $batchSize;

# if $cleanHeaders, don't print ##contig headers except for regular 
# chroms (\d+, X, Y, M/MT)
my $cleanHeaders = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = 'Merge several GVCFs into a single multi-sample GVCF, printed on STDOUT.
This software was developed and tested for merging single-sample GVCFs produced by 
Strelka (2.9.10), GATK (4.1.8.1), elPrep (5.0.2) or DeepVariant (1.4.0); and for 
multi-sample GVCFs produced by itself.
However it tries to be very defensive, so it should detect and report any problem it 
has with the GVCFs you provide. If this happens please report the issues so we can 
fix them.
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--filelist : file containing a list of GVCF filenames to merge, including paths, one per line
--tmpdir : subdir where tmp files will be created (on a RAMDISK if possible), must not pre-exist 
         and will be removed after execution
--jobs ['.$jobs.'] : number of parallel jobs=threads to run
--batchSize [adaptive] : size of each batch, lower decreases RAM requirements and performance,
        defaults to an optimized adaptive strategy, you shouldn\'t need to specify this except
        if you are running out of RAM (suggested reasonable values: between 5000 and 100000)
--cleanheaders : don\'t print ##contig headers except for chroms \d+ and X,Y,M/MT (possibly chr-prefixed)
--help : print this USAGE';

# construct string with full command-line for adding to headers, must be
# done before GetOptions
my $addToHeader = "$0 ".join(" ",@ARGV);
chomp($addToHeader);
$addToHeader .= " > ".`readlink -f /proc/$$/fd/1` ;
chomp($addToHeader);
$addToHeader .= " 2> ".`readlink -f /proc/$$/fd/2` ;
chomp($addToHeader);
$addToHeader = "##mergeGVCFs=<commandLine=\"$addToHeader\">\n";

GetOptions ("filelist=s" => \$fileList,
            "tmpdir=s" => \$tmpDir,
            "jobs=i" => \$jobs,
            "batchSize=i" => \$batchSize,
            "cleanheaders" => \$cleanHeaders,
            "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($fileList) ||
    die "E $0: you MUST provide a list of GVCFs to merge in a file, with --filelist.\n$USAGE\n";
(-f $fileList) ||
    die "E $0: provided --filelist $fileList doesn't exist or can't be read\n";

($tmpDir) || die "E $0: you must provide a tmpDir\n";
(-e $tmpDir) && 
    die "E $0: found argument $tmpDir but it already exists, remove it or choose another name.\n";
mkdir($tmpDir) || die "E $0: cannot mkdir tmpDir $tmpDir\n";

if ($jobs <= 2) {
    #  need one thread for eatTmpFiles and at least one worker thread
    warn "W $0: you set jobs=$jobs but we need at least 2 jobs, setting jobs=2 and proceeding\n";
    $jobs = 2;
}

if ($batchSize) {
    # disable adaptive strategy
    $batchSizeAdaptive = 0;
}
else {
    # default inital value
    $batchSize = 5000;
}


#############################################
## fill @infiles

# array of filehandles, open for reading each GVCF infile
my @infiles = ();

open(FILES, $fileList) ||
    die "E $0: cannot open provided argument $fileList\n";

while(my $file = <FILES>) {
    chomp($file);
    # skip blank lines
    ($file =~ /^\s*$/) && next;
    (-f $file) || die "E $0: infile $file from fileList $fileList doesn't exist\n";
    if ($file =~ /\.vcf$/) {
        open(my $infile, $file) || die "E $0: cannot open infile $file for reading\n";
        push(@infiles, $infile);
    }
    elsif ($file =~ /\.vcf\.gz$/) {
        # allow $jobs/8 threads for each infile, they should auto-regulate since
        # each bgzip process will have to wait for it's output to be eaten by the
        # pipe, and we stall worker threads when needed
        my $bgunzipJobs = int($jobs/8);
        ($bgunzipJobs > 0) || ($bgunzipJobs = 1);
        open(my $infile, "bgzip -c -d -\@$bgunzipJobs $file |") || die "E $0: cannot bgzip-open infile $file for reading\n";
        push(@infiles, $infile);
    }
    else {
        die "E $0: infile $file from fileList $fileList does not end in .vcf or .vcf.gz, why?\n";
    }
}

close(FILES);


#############################################
# deal with headers and fill @numSamples

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - STARTING TO WORK\n";

# same indexes as @infiles, value is the number of samples in the infile
my @numSamples;

# end of header == current command-line and #CHROM line . This allows to
# print the sample-column IDs from all files.
my $headerEnd = "$addToHeader";

# copy header from first file
my $infile = $infiles[0];
while(my $line = <$infile>) {
    if (($line =~ /^##INFO/) && ($line !~ /^##INFO=<ID=END,/) && ($line !~ /^##INFO=<ID=BLOCKAVG_min30p3a,/)) {
        next;
    }
    elsif (($cleanHeaders) && ($line =~ /^##contig=<ID=([^,]+),/)) {
        my $contig = $1;
        $contig =~ s/^chr//;
        if (($contig =~ /^\d+$/) || ($contig =~ /^[XYM]$/) || ($contig =~ /^MT$/)) {
            print $line;
        }
        # otherwise this is not a regular chrom, don't print
    }
    elsif ($line =~ /^##/) {
        print $line;
    }
    elsif ($line =~ /^#CHROM/) {
        # add full command-line to headers
        chomp($line);
        $headerEnd .= $line;
        # VCF has 9 columns in addition to the data columns
        my @f = split(/\t/,$line);
        push(@numSamples, scalar(@f) - 9);
        last;
    }
    else {
        die "E $0: parsing header from first infile of $ARGV[0], found bad line:\n$line";
    }
}

# skip headers from other infiles but grab ##mergeGVCFs and sample ids
foreach my $i (1..$#infiles) {
    my $infile = $infiles[$i];
    while(my $line = <$infile>) {
        if ($line =~ /^##mergeGVCFs=/) {
            print $line;
        }
        elsif ($line =~ /^##/) {
            #NOOP, skip
        }
        elsif ($line =~ /^#CHROM/) {
            # grab sample ids
            chomp($line);
            my @f = split(/\t/,$line);
            push(@numSamples, scalar(@f) - 9);
            ($line =~ /\tFORMAT(\t.+)$/) ||
                die "E $0: parsing CHROM header line from file $i:\n$line\n";
            $headerEnd .= $1;
            last;
        }
        else {
            die "E $0: parsing header from infile $ARGV[$i], found bad line:\n$line";
        }
    }
}

print "$headerEnd\n";

# sanity
(@numSamples == @infiles) || 
    die "E $0: numSamples and infiles disagree...\n@numSamples\n@infiles\n";

# flush stdout before starting our eatTmpFiles job
STDOUT->flush();

#############################################
# done with headers, now deal with bodies...

# create fork manager
my $pm = new Parallel::ForkManager($jobs);

# spawn a child process that waits for workers to finish producing batches,
# and prints the tmpfiles to stdout in correct order, cleaning up behind itself
if (! $pm->start) {
    &eatTmpFiles($tmpDir);
    $pm->finish;
}

# array of "lines", a line is an arrayref to a tab-splitted (up to DATA columns,
# ie all DATA columns are kept as a single un-split string as element [9]) line
# from infile, one line (max) per infile, each line has been parsed from the infile
#  but didn't belong to the current batch and will be dealt with in the next batch
my @startNextBatch;

# initialize with first non-filtered dataline from each file
foreach  my $i (0..$#infiles) {
    my $infile = $infiles[$i];
    my $lineR = &grabNextLine($infile);
    # every infile must have at least one non-filtered dataline
    ($lineR) || 
        die "E $0: infile $i doesn't have ANY non-filtered dataline?? Something must be wrong";
    push(@startNextBatch, $lineR);
}

# need $batchNum for the tmpFiles
my $batchNum = 0;

# index of first file that still has at least one line
my $firstFile = 0;

# key== a CHROM that we believe we are done with, value==1
# this allows to catch cases where the $firstFile doesn't contain
# any line for a chrom but has lines for subsequent chroms.
# Currently if this happens we just die, implementing a correct
# processing of this corner-case isn't trivial and I don't expect 
# it to happen...
my %chromsDone;

# timestamp for batchSizeAdaptive strategy
my $timestampAdaptive = time();

while ($firstFile <= $#infiles) {
    # precondition: we must have a line for file $firstFile in @startNextBatch
    $batchNum++;
    if (($batchSizeAdaptive) && ($batchNum % ($jobs-1) == 0) && ($batchNum >= $jobs)) {
        # jobs-1 because we want #workerThreads, excluding the eatTmpFiles job
        # $batchNum >= $jobs so we don't adjust on first batch of workers
        my $newTime = time();
        my $elapsed = $newTime - $timestampAdaptive;
        if ($elapsed < $batchTimeLow) {
            # increase batchSize by factor (1.2 * $btLow / $elapsed)
            $batchSize = int(1.2 * $batchSize * $batchTimeLow / $elapsed);
            # however, never exceed $maxIV, otherwise when looping on 1..$batchSize we
            # get an error "Range iterator outside integer range"
            my $maxIV = -1 >> 1;
            ($batchSize > $maxIV) && ($batchSize = $maxIV);
            $now = strftime("%F %T", localtime);
            warn "I $now: $0 - batchNum=$batchNum, adjusting batchSize up to $batchSize\n";
        }
        elsif ($elapsed > $batchTimeHigh) {
            # decrease batchSize by factor $btHigh / (1.2 * $elapsed) 
            $batchSize = int($batchSize * $batchTimeHigh / $elapsed / 1.2);
            $now = strftime("%F %T", localtime);
            warn "I $now: $0 - batchNum=$batchNum, adjusting batchSize down to $batchSize\n";
        }
        $timestampAdaptive = $newTime;
    }

    # chrom to deal with in the current batch, grab it from first file
    my $thisChr = $startNextBatch[$firstFile]->[0] ;
    ($chromsDone{$thisChr}) &&
        die "E $0: found line with chrom $thisChr in file $firstFile, but we thought we had processed all lines for this chrom. ".
        "Expected cause: some file < $firstFile has NO lines for a chrom (probably just before $thisChr). ".
        "Fix the code if you want to deal with this corner-case.";
    ($firstFile > 0) &&
        warn "W $0: surprisingly one or more infiles (indexes 0-".($firstFile-1).
        ") don't have any variants for chrom $thisChr, maybe check that they aren't truncated. Proceeding anyways.\n";

    # array of refs (one per infile) to arrays of (lines == arrayrefs)
    my @batchToMerge;
    # move stored line from first file to @batchToMerge
    $batchToMerge[$firstFile] = [$startNextBatch[$firstFile]];
    $startNextBatch[$firstFile] = undef;

    # smallest position to deal with in next batch, or 0 if all remaining
    # lines for current chrom must be eaten (ie next batch is another chrom or 
    # this is the last batch). Initialize with some non-zero value.
    my $posNextBatch = 1;

    # start with first file, so we can set $posNextBatch
    my $infile = $infiles[$firstFile];
    foreach my $j (1..$batchSize) {
        if (my $lineR = &grabNextLine($infile)) {
            my $chr = $lineR->[0];
            if ($chr eq $thisChr) {
                &addLineToBatch($batchToMerge[$firstFile], $lineR);
            }
            else {
                # $lineR is from another chrom, will be for next batch
                $posNextBatch = 0;
                $startNextBatch[$firstFile] = $lineR;
                $chromsDone{$thisChr} = 1;
                last;
            }
        }
        else {
            # no more non-filtered lines in $infile
            $posNextBatch = 0;
            last;
        }
    }
    # 3 possibilities:
    # - we ran out of lines -> $posNextBatch==0
    # - we read a line with a different chrom -> $posNextBatch==0 again
    # - otherwise we have one too many lines in @{$batchToMerge[$firstFile]},
    #   move it to @startNextBatch and set $posNextBatch
    if ($posNextBatch) {
        my $lineR = pop(@{$batchToMerge[$firstFile]});
        $posNextBatch = $lineR->[1];
        $startNextBatch[$firstFile] = $lineR;
    }

    # Now deal with all files except the first:
    # move relevant previously stored lines of non-first files to @batchToMerge
    foreach my $i ($firstFile+1..$#infiles) {
        if (my $lineR = $startNextBatch[$i]) {
            my $chr = $lineR->[0];
            my $pos = $lineR->[1];
            if (($chr eq $thisChr) && (($pos < $posNextBatch) || ($posNextBatch == 0))) {
                # stored line for file $i is from good chrom and pos, we will move it to
                # @batchToMerge, but if it's a long non-var block that extends beyond 
                # $posNextBatch we also have to continue the block for the next batch.
                # splitBlockLine does it all
                $startNextBatch[$i] = &splitBlockLine($lineR, $posNextBatch-1);
                $batchToMerge[$i] = [$lineR];
            }
            # else: file $i doesn't have any more lines for this chrom before $posNextBatch,
            # leave its next line in @startNextBatch ie NOOP
        }
    }

    # fill @batchToMerge for non-first files from the infiles themselves
    foreach my $i ($firstFile+1..$#infiles) {
        # if we already have a line in startNextBatch, this file doesn't
        # have any more lines for $thisChr before $posNextBatch
        ($startNextBatch[$i]) && next;
        # otherwise:
        my $infile = $infiles[$i];
        while (my $lineR = &grabNextLine($infile)) {
            my $chr = $lineR->[0];
            my $pos = $lineR->[1];
            if (($chr eq $thisChr) && (($pos < $posNextBatch) || ($posNextBatch == 0))) {
                $startNextBatch[$i] = &splitBlockLine($lineR, $posNextBatch-1);
                &addLineToBatch($batchToMerge[$i], $lineR);
                # if we produced a continuation blockline, this was the last line for infile $i
                ($startNextBatch[$i]) && last;
            }
            else {
                # wrong chrom or position too big, store for next batch
                $startNextBatch[$i] = $lineR;
                last; # last line for infile $i
            }
        }
    }

    # increment $firstFile if needed, for next batch
    while ((!$startNextBatch[$firstFile]) && ($firstFile <= $#infiles)) {
        $firstFile++;
    }

    # OK we can merge the batch and print the result to a tmpfile, but wait a bit if we already
    # have a lot of *done tmpfiles, otherwise in some scenarios we saturate tmpDir
    while(1) {
        my @doneFiles = glob("$tmpDir/*_done");
        (@doneFiles <= 1 + int($jobs / 4)) && last;
        sleep(10);
    }
    
    $pm->start && next;
    # if you change $tmpOut or $tmpOutFlag you MUST EDIT &eatTmpFiles()
    my $tmpOut = "$tmpDir/$batchNum.g.vcf";
    open(my $outFH, "> $tmpOut") || die "E $0: cannot open $tmpOut for writing\n";
    &mergeBatchOfLines(\@batchToMerge, \@numSamples, $outFH);
    close($outFH) || die "E $0: cannot close $tmpOut after writing ($!), is your tmpDir full?\n";
    # done building this tmpfile, create flag-file
    my $tmpOutFlag = $tmpOut."_done";
    open(OUTFLAG, "> $tmpOutFlag") || die "E $0: cannot open flagfile $tmpOutFlag for writing\n";
    print OUTFLAG "$batchNum\n";
    close(OUTFLAG) || die "E $0: cannot close $tmpOutFlag ($!), is your tmpDir full?\n";
    $pm->finish;
}

# some children are still processing batches, but we know the last
# batchNum that will ever exist, tell &eatTmpFiles() so it can exit
# (of course if you change $tmpOutLast you have to edit &eatTmpFiles)
my $tmpOutLast = "$tmpDir/lastBatch";
open(OUTLAST, "> $tmpOutLast") || die "E $0: cannot open tmp-last-file $tmpOutLast for writing\n";
print OUTLAST "$batchNum\n";
close(OUTLAST) || die "E $0: cannot close $tmpOutLast ($!), is your tmpDir full?\n";

$pm->wait_all_children;
# sanity: 
foreach my $i (0..$#infiles) {
    ($startNextBatch[$i]) && 
        die "E $0: All done but startNextBatch not empty for file $i:\n".join("\t",@{$startNextBatch[$i]})."\n";
}
foreach my $infile (@infiles) {
    close($infile);
}

$now = strftime("%F %T", localtime);
rmdir($tmpDir) || 
    die "E $now: $0 - all done but cannot rmdir tmpDir $tmpDir, why? $!\n";
warn "I $now: $0 - ALL DONE\n";



#############################################
## subs


###############
# grabNextLine:
# take a single arg: a filehandle open for reading (GVCF, after headers);
# return the next "line" that doesn't have a %filtersApplied in FILTER
# and that isn't a useless "GT:PGT:PID:PS" line from GATK,
# or undef if no such line exists.
# Also the variants are normalized: 
# trailing bases common to REF and all ALTs are removed, 
# leading bases common to REF and all ALTs are also removed and POS is 
# adjusted accordingly.
# A "line" is an arrayref, result of tab-splitting the actual line up to
# DATA columns, with FILTER cleared ('.') and INFO also cleared except if
# line is a non-var block (ie has END=).
# The last array element has all DATA columns in a single string.
sub grabNextLine {
    (@_ == 1) || die "E $0: grabNextLine needs 1 arg.\n";
    my ($infile) = @_;
  LINE:
    while (my $li = <$infile>) {
        chomp($li);
        # split into 10 fields, last one has all DATA columns
        my @line = split(/\t/, $li, 10);
        # ignore GATK "GT:PGT:PID:PS" lines, these are useless
        ($line[8] eq "GT:PGT:PID:PS") && (next LINE);
        # filter bad lines
        my @filters = split(/;/,$line[6]);
        foreach my $f (@filters) {
            (defined $filtersApplied{$f}) && (next LINE);
        }
        # if we get here no filters apply, clear FILTER
        $line[6] = ".";
        # also clear INFO if not in a non-var block
        ($line[7] =~ /^END=/) || ($line[7] = ".");

        # normalize variants:
        # most REFs are a single base => cannot be normalized -> test first
        if (length($line[3]) >= 2) {
            my $ref = $line[3];
            my @alts = split(/,/,$line[4]);
            # never normalize <NON_REF>, * or <*> (the DV equivalent of gatk's <NON_REF>) : if they are here,
            # store their indexes in @alts and splice them out (trick: start from the end)
            my ($nonrefi,$stari,$dvStari) = (-1,-1,-1);
            foreach my $i (reverse(0..$#alts)) {
                if ($alts[$i] eq '<NON_REF>') {
                    $nonrefi = $i;
                    splice(@alts,$i,1);
                }
                elsif ($alts[$i] eq '<*>') {
                    $dvStari = $i;
                    splice(@alts,$i,1);
                }
                elsif ($alts[$i] eq '*') {
                    $stari = $i;
                    splice(@alts,$i,1);
                }
            }

            # 1. if length >= 2 for REF and all ALTS, and if REF and all ALTs have 
            #    common ending bases, remove them (keeping at least 1 base everywhere).
            while ($ref =~ /\w(\w)$/) {
                # ref has at least 2 chars
                my $lastRef = $1;
                my $removeLast = 1;
                foreach my $alt (@alts) {
                    if ($alt !~ /\w$lastRef$/) {
                        # this alt is length one or doesn't end with $lastRef
                        $removeLast = 0;
                        last;
                    }
                }
                if ($removeLast) {
                    # OK remove last base from REF and all @alts
                    ($ref =~ s/$lastRef$//) || 
                        die "E $0: WTF can't remove $lastRef from end of $ref\n";
                    foreach my $i (0..$#alts) {
                        ($alts[$i] =~ s/$lastRef$//) || 
                            die "E $0: WTF can't remove $lastRef from end of alt $i == $alts[$i]\n";
                    }
                }
                else {
                    # can't remove $lastRef, get out of while loop
                    last;
                }
            }

            # 2. if length >= 2 for REF and all ALTS, and if REF and all ALTs have 
            #    common starting bases, remove them (keeping at least 1 base everywhere)
            #    and adjust POS.
            while ($ref =~ /^(\w)\w/) {
                my $firstRef = $1;
                my $removeFirst = 1;
                foreach my $alt (@alts) {
                    if ($alt !~ /^$firstRef\w/) {
                        $removeFirst = 0;
                        last;
                    }
                }
                if ($removeFirst) {
                    ($ref =~ s/^$firstRef//) || 
                        die "E $0: WTF can't remove $firstRef from start of $ref\n";
                    foreach my $i (0..$#alts) {
                        ($alts[$i] =~ s/^$firstRef//) || 
                            die "E $0: WTF can't remove $firstRef from start of alt $i == $alts[$i]\n";
                    }
                    $line[1]++;
                }
                else {
                    last;
                }
            }

            # place <NON_REF> and/or * or <*> back where they belong: need to splice
            # in correct order, smallest index first
            # NOTE: dvStari (deepVariant only) and stari/nonrefi (GATK/elPrep) are exclusive
            if ($stari != -1) {
                if ($nonrefi == -1) {
                    splice(@alts, $stari, 0, '*');
                }
                elsif ($nonrefi > $stari) {
                    splice(@alts, $stari, 0, '*');
                    splice(@alts, $nonrefi, 0, '<NON_REF>');
                }
                else {
                    # nonref needs to be spliced back in first, then *
                    splice(@alts, $nonrefi, 0, '<NON_REF>');
                    splice(@alts, $stari, 0, '*');
                }
            }
            elsif ($nonrefi != -1) {
                splice(@alts, $nonrefi, 0, '<NON_REF>');
            }
            elsif ($dvStari != -1) {
                splice(@alts, $dvStari, 0, '<*>');
            }

            $line[3] = $ref;
            $line[4] = join(',',@alts);
        }
        
        return(\@line);
  }
    # no more good lines in infile, return undef
    return();
}

###############
# addLineToBatch: args are:
# - $batchR, a ref to an array of arraylines;
# - $lineR, a ref to an arrayline.
# This function "pushes" $lineR onto @$batchR, but does the right thing
# when $lineR->POS isn't simply greater than the POS of the line at the
# top of $batchR.
# This happens when &grabNextLine shifted POS of the previous line,
# or also due to some Strelka bugs (sometimes a non-var call is made
# (whether individually or in a block) and then in the next line an indel
# call is made with the same POS)
# Precondition, caller must ensure it (not checked here!):
# $lineR->chrom is the same as chrom in batchR
sub addLineToBatch {
    (@_ == 2) || die "E $0: addLineToBatch needs 2 args.\n";
    my ($batchR,$lineR) = @_;

    # if batchR is empty, fast return (can happen due to recursion below)
    if (! @$batchR) {
        push(@$batchR, $lineR);
        return();
    }

    my $prevLineR = $batchR->[$#$batchR];
    
    if ($lineR->[1] > $prevLineR->[1]) {
        # newPOS > prevPOS, just make prevLine end at newPos-1 if it was a block
        &splitBlockLine($prevLineR, $lineR->[1] - 1);
        push(@$batchR, $lineR);
    }
    elsif ($lineR->[1] == $prevLineR->[1]) {
        if (($lineR->[4] eq '.') || ($lineR->[4] eq '<NON_REF>') || ($lineR->[4] eq '<*>')) {
            # if lineR is non-var ignore it at POS, but make it start at POS+1
            # if it's a block (just skip it if it's not a block)
            $lineR = &splitBlockLine($lineR,$lineR->[1]);
            ($lineR) && (push(@$batchR,$lineR));
        }
        elsif(($prevLineR->[4] eq '.') || ($prevLineR->[4] eq '<NON_REF>') || ($prevLineR->[4] eq '<*>')) {
            # prevLineR was non-var, replace it with lineR
            pop(@$batchR);
            push(@$batchR,$lineR);
        }
        else {
            # this shouldn't happen, it would mean infile has contradictory 
            # calls at the same POS...
            # we could try to deal with it (it DOES happen with strelka...),
            # but it's rare and usually corresponds to 2 HET calls on different
            # ALT alleles (ie would be 1/2 calls if we merged them). These calls
            # are typically ignored downstream (they go in the OTHER column).
            # So, for now just keep the first line, ie NOOP

            # following warn was used to analyze the problem, don't use it in production
            #warn "W $0: addLineToBatch, 2 lines with ALTs and same POS, ignoring the second line:\n".
            # join("\t",@$prevLineR)."\n".join("\t",@$lineR)."\n\n";
        }
    }
    else {
        # newPOS < prevPOS, we can switch the 2 lines and recurse,
        # but how far back do we have to go?
        # easy to code with recursive function, but I sure hope
        # this doesn't happen too often!
        warn "I $0: addLineToBatch is recursing, around $lineR->[0] $lineR->[1]\n";
        pop(@$batchR);
        &addLineToBatch($batchR, $lineR);
        &addLineToBatch($batchR, $prevLineR);
    }
}


###############
# splitBlockLine: args are:
# - $lineR, a ref to an arrayline;
# - $newEnd, an int which must be >= $lineR->POS or < 0.
# -> if $lineR doesn't contain a non-variant block: NOOP & return undef
# -> elsif $newEnd < 0: NOOP & return undef
# -> elsif $newEnd >= $lineR->END: NOOP & return undef
# -> else return a ref to a freshly created arrayline, with same 
#    content as $lineR except POS becomes $newEnd+1 and REF becomes 'N';
#    and also modify $lineR: END becomes $newEnd.
sub splitBlockLine {
    (@_ == 2) || die "E $0: splitBlockLine needs 2 args.\n";
    my ($lineR,$newEnd) = @_;

    # if $newEnd < 0, NOOP
    ($newEnd >= 0) || (return());
    # if $lineR isn't a block, NOOP
    ($lineR->[7] =~ /^END=(\d+)/) || (return());
    my $oldEnd = $1;
    # if $lineR is a block but ends before or at $newEnd, NOOP
    ($newEnd < $oldEnd) ||  (return());

    # sanity: can't have newEnd positive but < POS
    ($newEnd < $lineR->[1]) &&
        die "E $0: splitBlockLine, newEnd $newEnd smaller than POS in:\n".join("\t",@$lineR)."\n";

    # else make fresh line, initially copy content of $lineR
    my @newLine = @$lineR ;
    # replace POS with $newEnd+1 and REF with N
    $newLine[1] = $newEnd + 1;
    $newLine[3] = 'N';
    # modify END in $lineR
    ($lineR->[7] =~ s/^END=$oldEnd/END=$newEnd/) ||
        die "E $0: splitBlockLine, cannot subst END $oldEnd with $newEnd in:\n".join("\t",@$lineR)."\n";
    return(\@newLine);
}


###############
# mergeBatchOfLines: args are
# $batchToMergeR == arrayref of refs (one per infile) to arrays of 
#    "lines", each line is actually an arrayref to the tab-splitted line
# ref to @numSamples == same indexes as @infiles, value is the number
#    of samples in the infile
# $outFH, a filehandle open for writing
#
# Preconditions:
# lines from the first infile define the start and end chrom:pos for 
# this batch; 
# -> All lines in the batch have the same chrom
# -> All lines from any infile with pos >= the first file's starting
#   pos in this batch, and with chrom:pos < the first file's starting
#   chrom:pos in the NEXT batch, must be in this batch
# -> Any line storing a non-variant block has END < the first file's
#   starting chrom:pos in the NEXT batch (splitBlockLine facilitates this)
sub mergeBatchOfLines {
    (@_ == 3) || die "E $0: mergeBatchOfLines needs 3 args.\n";
    my ($batchToMergeR, $numSamplesR, $outFH) = @_;

    # number of files that still have at least one line to merge
    my $filesNotFinished = 0;
    foreach my $i (0..$#infiles) {
        ($batchToMergeR->[$i]->[0]) && ($filesNotFinished++);
    }

    while ($filesNotFinished) {
        # @nextToMerge: "array-lines" for the lines that have
        # the "smallest" position, ie that must be merged next.
        # same indexes as @$batchToMergeR.
        my @nextToMerge ;
        # current smallest pos
        my $smallestPos = 0;
        # longest REF among the toMerge lines
        my $longestRef = "";
        # FORMAT keys that appear in the toMerge lines
        my %longestFormat = ();
        # FORMAT strings that appear in toMerge lines
        my %seenFormats = ();
        # will hold the smallest END if lines to merge are all non-variant blocks, POS otherwise
        my $nonVarBlockEnd;

        foreach my $i (0..$#infiles) {
            # ignore files without remaining lines in this batch
            ($batchToMergeR->[$i]->[0]) || next;
            
            # $fieldsR just to simplify code
            my $fieldsR = $batchToMergeR->[$i]->[0];
            my ($pos,$ref,$format) = @$fieldsR[1,3,8];
            if (($pos < $smallestPos) || (! $smallestPos)) {
                # this is the new smallest
                @nextToMerge = ();
                $nextToMerge[$i] = $fieldsR;
                # update $smallestPos at the end, we need previous value for END=
                $longestRef = $ref;
                %longestFormat = ();
                foreach my $fkey (split(/:/,$format)) {
                    $longestFormat{$fkey} = 1;
                }
                %seenFormats = ();
                $seenFormats{$format} = 1;
                if ($fieldsR->[7] =~ /^END=(\d+)/) {
                    my $newBlockEnd = $1;
                    # examine the previous $smallestPos here!
                    if ((! $smallestPos) || ($newBlockEnd < $smallestPos)) {
                        # there was no previous smallestPos, or it was beyond
                        # the end of this block -> we can use the whole block
                        $nonVarBlockEnd = $newBlockEnd ;
                    }
                    else {
                        # we had a position starting later but within this block,
                        # we have to end this block just before
                        $nonVarBlockEnd = $smallestPos - 1;
                    }
                }
                else {
                    # this is not a blockline
                    $nonVarBlockEnd = $pos;
                }
                # in any case update $smallestPos
                $smallestPos = $pos;
            }
            elsif ($pos == $smallestPos) {
                $nextToMerge[$i] = $fieldsR;
                # replace longestRef if new is longer or if old was 'N'
                if ((length($ref) > length($longestRef)) || ($longestRef eq 'N')) {
                    $longestRef = $ref;
                }
                if (! $seenFormats{$format}) {
                    foreach my $fkey (split(/:/,$format)) {
                        $longestFormat{$fkey} = 1;
                    }
                    $seenFormats{$format} = 1;
                }
                if ($fieldsR->[7] =~ /^END=(\d+)/) {
                    my $thisEnd = $1;
                    ($thisEnd < $nonVarBlockEnd) && ($nonVarBlockEnd = $thisEnd);
                }
                else {
                    $nonVarBlockEnd = $pos;
                }
            }
            # remaining cases are all $pos > $smallestPos, but...
            elsif ($pos <= $nonVarBlockEnd) {
                # current file has "larger" position but selected is in a non-var block, we must stop
                # the non-var block just before current file's position
                $nonVarBlockEnd = $pos - 1;
            }
        }

        # remove the chosen lines from $batchToMergeR and decrement 
        # $filesNotFinished if we finished a file, except for non-variant
        # block lines where we only remove the lines if we eat the block
        # to its END (otherwise keep line but fix POS and REF, and also fix
        # END in the nextToMerge copy)
        foreach my $i (0..$#nextToMerge) {
            if ($nextToMerge[$i]) {
                my $continueBlockR = &splitBlockLine($nextToMerge[$i], $nonVarBlockEnd);
                if ($continueBlockR) {
                    $batchToMergeR->[$i]->[0] = $continueBlockR;
                }
                else {
                    shift(@{$batchToMergeR->[$i]});
                    ($batchToMergeR->[$i]->[0]) || ($filesNotFinished--);
                }
            }
        }

        if ($smallestPos < $nonVarBlockEnd) {
            # if we get here we are in a non-var block in EVERY file, otherwise
            # $nonVarBlockEnd would have been set to $pos
            # make sure all blocks have the same FORMAT
            (keys(%seenFormats) == 1) ||
                die "E $0: in mergeBatchOfLines I want to merge non-var blocks at pos $smallestPos but I have several formats:".
                keys(%seenFormats);
            print($outFH &mergeLinesNonVarBlock(\@nextToMerge,$numSamplesR, $longestRef));
        }
        else {
            # build array from %longestFormat, respecting the order specified 
            # in @maxFormatSorted, and ignoring MIN_DP since we are not in a non-var block
            my @maxFormatSorted = ('GT','AF','DP','FT','GQ','GQX','DPF','DPI','AD','ADF','ADR','VAF','SB','PL','PS','PGT','PID');
            my @longestFormat = ();
            foreach my $fkey (@maxFormatSorted) {
                ($longestFormat{$fkey}) && (push(@longestFormat,$fkey));
            }
            print($outFH &mergeLines(\@nextToMerge, $numSamplesR, $longestRef, \@longestFormat));
        }
    }
}


###############
# mergeLines: merge (up to) one line per infile, args are:
# arrayref of refs to array-lines, one array-line per infile, 
#    all array-lines have same chrom:pos or are undef
# ref to @numSamples == same indexes as @infiles, value is the number of samples in the infile
# $longestRef wich contains the longest REF among the strings
# ref to @longestFormat: array of all FORMAT keys present in at least one line 
#    (ignoring MIN_DP but adding FT if it's needed)
# Return the string to print, result of merging the lines.
# NOTE: arrays referenced in $toMergeR WILL BE modified (ALTs are changed)
# Preconditions: 
# - MIN_DP must be AFTER DP in FORMAT if it's there (this is checked)
# - any non-variant blockline should have END==POS (but we don't actually
#   examine or use END, so this doesn't matter)
sub mergeLines {
    (@_ == 4) || die "E $0: mergeLines needs 4 args.\n";
    my ($toMergeR, $numSamplesR, $longestRef, $longestFormatR) = @_;

    ####################################
    # STEP 1: fix ALTs in @$toMergeR so they correspond to a REF == $longestRef,
    # and fill %newAlts: 
    # key == a new ALT as it appears in the fixed @$toMergeR, and as it will appear in 
    # the returned merged string;
    # value will be 1 initially, but becomes the alt's index in @newAlts during STEP 2
    my %newAlts;
    foreach my $fileIndex (0..$#$toMergeR) {
        ($toMergeR->[$fileIndex]) || next;
        my ($ref,$alts) = @{$toMergeR->[$fileIndex]}[3,4];
        # if we're in a non-variant position in this file: $ref might be 'N' and anyways
        # no extension of alts is needed, but <NON_REF> (GATK) or <*> (DV) must still go in newAlts
        ($alts eq '.') && next;
        if (($alts eq '<NON_REF>') || ($alts eq '<*>')) {
            $newAlts{$alts} = 1;
            next;
        }
        if ($longestRef ne $ref) {
            ($longestRef =~ /^$ref(\w+)$/) || 
                die "E $0: longestRef $longestRef doesn't start with ref $ref (file $fileIndex), impossible\n".
                join("\t",@{$toMergeR->[$fileIndex]})."\n";
            my $extraBases = $1;
            my @fixedAlts = ();
            foreach my $thisAlt (split(/,/,$alts)) {
                my $fixedAlt = $thisAlt;
                # never extend <NON_REF> or * or <*>
                if (($fixedAlt ne '<NON_REF>') && ($fixedAlt ne '*') && ($fixedAlt ne '<*>')) {
                    $fixedAlt .= $extraBases;
                }
                $newAlts{$fixedAlt} = 1;
                push(@fixedAlts, $fixedAlt);
            }
            # now fix the actual line: replace ALTs, preserve order
            $toMergeR->[$fileIndex]->[4] = join(',',@fixedAlts);
        }
        else {
            # else $ref is the longestRef, just update %newAlts
            foreach my $thisAlt (split(/,/,$alts)) {
                $newAlts{$thisAlt} = 1;
            }
        }
    }

    ####################################
    # STEP 2 (NO LOOP): build @newAlts, the list of fixed ALTs sorted appropriately:
    # by increasing length, and at equal length in alphabetical order (so we are
    # deterministic), always keeping '*' and <NON_REF> and <*> last (in that order) if present
    my ($starPresent,$nonrefPresent,$dvStarPresent) = (0,0,0);
    ($newAlts{'*'}) && ($starPresent = 1) && (delete($newAlts{'*'}));
    ($newAlts{'<NON_REF>'}) && ($nonrefPresent = 1) && (delete($newAlts{'<NON_REF>'}));
    ($newAlts{'<*>'}) && ($dvStarPresent = 1) && (delete($newAlts{'<*>'}));

    my @newAlts = sort {(length($a) <=> length($b)) || ($a cmp $b)} keys(%newAlts);
    # add back *, NON_REF and <*> if they were here
    ($starPresent) && (push(@newAlts, '*'));
    ($nonrefPresent) && (push(@newAlts, '<NON_REF>'));
    ($dvStarPresent) && (push(@newAlts, '<*>'));
    # if no ALTs, use '.' as a bogus ALT
    (@newAlts) || push(@newAlts, '.');

    # also update %newAlts: key==ALT, value becomes the index of that alt in @newAlts
    foreach my $i (0..$#newAlts) {
        $newAlts{$newAlts[$i]} = $i;
    }

    ####################################
    # STEP 3 (NO LOOP): start building line that will be returned
    my $toPrint;
    # grab CHROM POS from first non-null file
    my $firstNonNull = 0;
    while(! $toMergeR->[$firstNonNull]) {
        $firstNonNull++;
    }
    $toPrint = $toMergeR->[$firstNonNull]->[0]."\t".$toMergeR->[$firstNonNull]->[1];
    # ID REF
    $toPrint .= "\t.\t".$longestRef;
    # ALT
    $toPrint .= "\t".join(',',@newAlts);
    # QUAL FILTER INFO
    $toPrint .= "\t.\t.\t.";
    # FORMAT: @$longestFormatR must NOT contain MIN_DP (checked below)
    $toPrint .= "\t".join(':',@$longestFormatR);
    # also need %formatIndex, key==FORMAT key, value==index of this key in @$longestFormatR
    my %formatIndex = ();
    foreach my $i (0..$#$longestFormatR) {
        ($longestFormatR->[$i] eq "MIN_DP") &&
            die "E $0: in mergeLines: longestFormat contains MIN_DP, it MUST NOT!\n";
        $formatIndex{$longestFormatR->[$i]} = $i;
    }

    ####################################
    # STEP 4: fix the data columns if needed and add them to $toPrint
    foreach my $fileIndex (0..$#$numSamplesR) {
        if (! $toMergeR->[$fileIndex]) {
            # no line, use blank data for all samples from this file
            foreach my $j (1..$numSamplesR->[$fileIndex]) {
                $toPrint .= "\t.";
            }
        }
        elsif (($toMergeR->[$fileIndex]->[3] eq $longestRef) && 
               ($toMergeR->[$fileIndex]->[4] eq join(',',@newAlts)) && 
               ($toMergeR->[$fileIndex]->[8] eq join(':',@$longestFormatR))) {
            # REF ALT FORMAT didn't change for this file, just copy the data
            $toPrint .= "\t".$toMergeR->[$fileIndex]->[9];
        }
        else {
            # no shortcut, have to examine and fix everything
            my @alts = split(/,/,$toMergeR->[$fileIndex]->[4]);
            # @altsNew2Old: same indexes as @newAlts,
            # value is the index of $newAlts[$i] in @alts, or undef if it's not there
            my @altsNew2Old;
            foreach my $i (0..$#alts) {
                # ignore '.', it cannot be referenced in @data
                ($alts[$i] eq '.') && next;
                $altsNew2Old[$newAlts{$alts[$i]}] = $i;
            }

            # make sure MIN_DP comes after DP if present
            ($toMergeR->[$fileIndex]->[8] =~ /:MIN_DP.*:DP/) &&
                die "E $0: in mergeLines MIN_DP is in a FORMAT but before DP! in file $fileIndex:\n".
                join("\t",@{$toMergeR->[$fileIndex]})."\n";

            my @format = split(/:/,$toMergeR->[$fileIndex]->[8]);

            # deal with each DATA column
            my @dataCols = split("\t",$toMergeR->[$fileIndex]->[9]);
            foreach my $j (0..$#dataCols) {
                my @data = split(/:/, $dataCols[$j]);

                # fixed DATA for this sample, one value per longestFormat key
                my @fixedData;

                foreach my $fi (0..$#format) {
                    #     GT gets adjusted values from %altsNew2Old
                    #     AF GQ GQX DPF DPI SB PS PGT PID FT don't change
                    #     DP gets MIN_DP value if it exists, otherwise doesn't change
                    #     AD ADF ADR VAF get 0 for new ALTs
                    #     PL gets correct new values (see explanation in the code)
                    if (!defined $data[$fi]) {
                        # data columns can be shorter than FORMAT when last values are unknown
                        last;
                    }
                    elsif ($format[$fi] eq "GT") {
                        if (($data[$fi] eq '.') || ($data[$fi] eq './.') || ($data[$fi] eq '.|.')) {
                            # nocall, just copy
                            $fixedData[$formatIndex{$format[$fi]}] = $data[$fi];
                        }
                        elsif ($data[$fi] =~ m~^(\d+)([/\|])(\d+)$~) {
                            # diploid, regular or phased
                            my ($geno1,$sep,$geno2) = ($1,$2,$3);
                            if ($geno1 != 0) {
                                $geno1 = 1 + $newAlts{$alts[$geno1-1]};
                            }
                            if ($geno2 != 0) {
                                $geno2 = 1 + $newAlts{$alts[$geno2-1]};
                            }
                            $fixedData[$formatIndex{$format[$fi]}] = "$geno1$sep$geno2";
                        }
                        elsif ($data[$fi] =~ m~^(\d+)$~) {
                            # hemizygous
                            my ($geno1) = ($1);
                            if ($geno1 != 0) {
                                $geno1 = 1 + $newAlts{$alts[$geno1-1]};
                            }
                            $fixedData[$formatIndex{$format[$fi]}] = "$geno1";
                        }
                        else {
                            die "E $0: trying to adjust GT $data[$fi] but can't parse it, in file $fileIndex, line:\n".
                                join("\t",@{$toMergeR->[$fileIndex]})."\n";
                        }
                    }

                    elsif (($format[$fi] eq "AF") || ($format[$fi] eq "GQ") || ($format[$fi] eq "GQX") || 
                           ($format[$fi] eq "DP") || ($format[$fi] eq "DPF") || ($format[$fi] eq "DPI") || 
                           ($format[$fi] eq "SB") || ($format[$fi] eq "FT") ||
                           ($format[$fi] eq "PS") || ($format[$fi] eq "PGT") || ($format[$fi] eq "PID") ) {
                        # just copy values
                        $fixedData[$formatIndex{$format[$fi]}] = $data[$fi];
                    }

                    elsif ($format[$fi] eq "MIN_DP") {
                        # we verified that if MIN_DP exists it is after DP in FORMAT,
                        # therefore if we find MIN_DP we can just use it for DP and 
                        # it will squash any pre-exising value
                        (defined $formatIndex{"DP"}) && ($fixedData[$formatIndex{"DP"}] = $data[$fi]);
                    }

                    elsif (($format[$fi] eq "AD") || ($format[$fi] eq "ADF") || ($format[$fi] eq "ADR") || ($format[$fi] eq "VAF")) {
                        my @ADs = split(/,/, $data[$fi]);
                        # copy REF AD and remove from @ADs (except for VAF where there's no REF value)
                        ($format[$fi] ne "VAF") && ($fixedData[$formatIndex{$format[$fi]}] = shift(@ADs).",");
                        # append ADs for all ALTs, inserting '0' except for 
                        # alts also present in currentLine
                        foreach my $newAltIndex (0..$#newAlts) {
                            if (defined $altsNew2Old[$newAltIndex]) {
                                $fixedData[$formatIndex{$format[$fi]}] .= $ADs[$altsNew2Old[$newAltIndex]].","; 
                            }
                            else {
                                # this newAlt is not in current line, use zero
                                $fixedData[$formatIndex{$format[$fi]}] .= "0,";
                            }
                        }
                        # remove trailing ','
                        chop($fixedData[$formatIndex{$format[$fi]}]);
                    }

                    elsif ($format[$fi] eq "PL") {
                        # ALTs that didn't exist got 0 counts, and the corresponding PLs 
                        # will get 255 values
                        # The VCF spec says the PL for x/y genotype is at index x + y*(y+1)/2
                        # Note that x and y here start at 1 for ALTs, 0 is the REF.
                        my @PLs = split(/,/, $data[$fi]);

                        my @newPLs;

                        # $badPLs for Strelka bug where sometimes PL doesn't have correct number of values
                        my $badPLs = 0;

                        # $x and $y start at -1 because we want them to be indexes of @altsNew2Old,
                        # using -1 for REF
                        foreach my $x (-1..$#newAlts) {
                            foreach my $y ($x..$#newAlts) {
                                # index in @newPLs where PL for X/Y will go,
                                # +1 to each because ALTs start at 1 in the VCF formula
                                my $plDest = ($x+1) + ($y+1)*($y+2)/2 ;
                                # index in @PLs where we grab the PL for X/Y
                                # careful: x and y may be in a different order in the source,
                                # in which case we must switch them
                                my $plSource = 0;
                                # $oldX and $oldY use the VCF convention => 0 is REF, ALTs start at 1
                                my ($oldX,$oldY) = (0,0);

                                if ($x == -1) {
                                    # X is REF -> NOOP, leave $oldX==0
                                }
                                elsif (defined $altsNew2Old[$x]) {
                                    $oldX = $altsNew2Old[$x] + 1;
                                }
                                else {
                                    # X unknown, use 255 as PL
                                    $newPLs[$plDest] = 255;
                                    next;
                                }

                                if ($y == -1) {
                                    # Y is REF -> NOOP
                                }
                                elsif (defined $altsNew2Old[$y]) {
                                    $oldY = $altsNew2Old[$y] + 1;
                                }
                                else {
                                    # Y unknown
                                    $newPLs[$plDest] = 255;
                                    next;
                                }

                                # switch if order is wrong
                                if ($oldX > $oldY) {
                                    my $tmp = $oldX;
                                    $oldX = $oldY;
                                    $oldY = $tmp;
                                }
                                
                                # calculate $plSource using VCF formula
                                $plSource = $oldX + $oldY*($oldY+1)/2;
                                # finally, fill newPLs
                                # $PLs[$plSource] should always exist according to spec, but Strelka has
                                # a bug where sometimes values are missing, and we can't know which 
                                # genotype the existing PLs were for...
                                # IOW in these cases the PL string really can't be interpreted!
                                # so we flag it here and then replace the whole PL value with comma-separated '.'
                                if (defined $PLs[$plSource]) {
                                    $newPLs[$plDest] = $PLs[$plSource];
                                }
                                else {
                                    $newPLs[$plDest] = '.';
                                    $badPLs = 1;
                                }
                            }
                        }
                        if ($badPLs) {
                            # wrong number of PL values in Strelka file, replace all values with '.'
                            $fixedData[$formatIndex{"PL"}] = '.';
                            $fixedData[$formatIndex{"PL"}] .= ",." x $#newPLs ; # we put one '.' already, need $# more
                        }
                        else {
                            $fixedData[$formatIndex{"PL"}] = join(',', @newPLs);
                        }
                    }

                    else {
                        die "E $0: unknown format key $format[$fi] found, add some code to deal with it!\n".
                            join("\t",@{$toMergeR->[$fileIndex]})."\n";
                    }
                }

                # now prune trailing dummy fields that may have been copied from infiles (but never GT)
                foreach my $i (reverse(1..$#fixedData)) {
                    if ($fixedData[$i] =~ /^[\.,]+$/) {
                        pop(@fixedData);
                    }
                    else {
                        last;
                    }
                }

                # if AD* / VAF is needed but wasn't found, use correct number of dummy '.'
                foreach my $fkey ("AD","ADF","ADR","VAF") {
                    if ((defined $formatIndex{$fkey}) && ($#fixedData >= $formatIndex{$fkey}) && (! $fixedData[$formatIndex{$fkey}])) {
                        # need one value for REF (except for VAF) + one value per newAlts allele, use dummy '.'
                        ($fkey ne 'VAF') && ($fixedData[$formatIndex{$fkey}] = '.,');
                        $fixedData[$formatIndex{$fkey}] .= ".," x scalar(@newAlts);
                        # remove last ','
                        chop($fixedData[$formatIndex{$fkey}]);
                    }
                }

                # if PL was asked for but not found: I don't want to pretend having 
                # PLs if I don't, don't want to bother coding it, and anyways PL
                # comes close to the end (so will rarely be needed & missing)
                # => pretend PL is single-valued

                # all remaining fields are single-valued, use dummy '.' if needed and missing
                foreach my $i (0..$#fixedData) {
                    (defined $fixedData[$i]) || ($fixedData[$i] = '.');
                }

                $toPrint .= "\t".join(':',@fixedData);
            }
        }
    }

    ####################################
    # STEP 5: all done! return result
    return("$toPrint\n");
}



###############
# special case of mergeLines: every line to merge is a non-variant block,
# and they all have the same POS, FORMAT and END (this is checked).
# -> return a new \n-terminated blockline with all samples merged
sub mergeLinesNonVarBlock {
    (@_ == 3) || die "E $0: mergeLinesNonVarBlock needs 3 args.\n";
    my ($toMergeR,$numSamplesR, $longestRef) = @_;

    # find first non-null line
    my $firstNonNull = 0;
    while(! $toMergeR->[$firstNonNull]) {
        $firstNonNull++;
    }

    # make a new non-variant block: grab start of line from first non-null file,
    # except for REF: use $longestRef because if any infile has non-N it will be there 
    my $toPrint = join("\t", @{$toMergeR->[$firstNonNull]}[0..2]);
    $toPrint .= "\t$longestRef";
    # grab ALT, QUAL, FILTER ('.'), INFO and FORMAT from first non-null file, all
    # files should have the same (this is checked for INFO and FORMAT)
    $toPrint .= "\t".join("\t", @{$toMergeR->[$firstNonNull]}[4..8]);

    # add data columns (even for files that don't have a line at this pos)
    foreach my $fileIndex (0..$#$numSamplesR) {
        if ($toMergeR->[$fileIndex]) {
            # check INFO
            ($toMergeR->[$fileIndex]->[7] eq $toMergeR->[$firstNonNull]->[7]) ||
                die "E $0: in mergeLinesNonVarBlock, INFO mismatch when toPrint=$toPrint\n";
            # check FORMAT
            ($toMergeR->[$fileIndex]->[8] eq $toMergeR->[$firstNonNull]->[8]) || 
                die "E $0: in mergeLinesNonVarBlock, FORMAT mismatch when toPrint=$toPrint\n";
            # print data columns (they are all in ->[9], we never split them)
            $toPrint .= "\t".$toMergeR->[$fileIndex]->[9];
        }
        else {
            # file $fileIndex doesn't have a line for this position, print 
            # "blank" data for each sample from infile
            foreach my $j (1..$numSamplesR->[$fileIndex]) {
                $toPrint .= "\t.";
            }
        }
    }
    $toPrint .= "\n";
    return($toPrint);
}


###############
# this function waits for flagfiles to be created and "eats" the 
# corresponding tmpFiles in order, starting at 1.
# "eating" means print to stdout and remove.
# we also watch for $tmpOutLast, a file that will tell us
# the last batch number to wait for
sub eatTmpFiles {
    (@_ == 1) || die "E $0: eatTmpFiles needs 1 arg.\n";
    my ($tmpD) = @_;

    # when created, $tmpOutLast will contain the number of the last batch
    my $tmpOutLast = "$tmpD/lastBatch";
    # $lastBatch remains undefined until we read it from $tmpOutLast
    my $lastBatch;

    # next batch to eat
    my $nextBatch = 1;

    while(1) {
        my $tmpOut = "$tmpD/$nextBatch.g.vcf";
        my $tmpOutFlag = $tmpOut."_done";

        if ((-e $tmpOutFlag)  && (! -z $tmpOutFlag)) {
            # next batch tmp file is ready
            open (IN, $tmpOut) || 
                die "E $0: in eatTmpFiles, flagfile $tmpOutFlag exists but cant open tmpout $tmpOut: $!\n";
            while(<IN>) {
                print;
            }
            close(IN);
            (unlink($tmpOut,$tmpOutFlag) == 2) ||
                die "E $0: in eatTmpFiles, done with files for batch $nextBatch but cannot unlink (both of) them: $!\n";

            if ($nextBatch % 10 == 0) {
                # log progress every 10 batches
                my $now = strftime("%F %T", localtime);
                warn "I $now: $0 - done printing results from batch $nextBatch\n";
            }
            $nextBatch++;
            next;
        }

        elsif ((-e $tmpOutLast) && (! -z $tmpOutLast)) {
            open (IN, $tmpOutLast) || 
                die "E $0: in eatTmpFiles, cannot open tmpOutLast $tmpOutLast although it exists: $!\n";
            $lastBatch = <IN>;
            chomp($lastBatch);
            close(IN);
            unlink($tmpOutLast) || 
                die "E $0: in eatTmpFiles, cannot unlink tmpOutLast $tmpOutLast: $!\n";
            next;
        }

        elsif ((defined $lastBatch) && ($lastBatch < $nextBatch)) {
            # all done, return so this process can finish
            return();
        }

        else {
            # wait a few seconds before looking again
            sleep(10);
            next;
        }
    }
}

