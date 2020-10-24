#!/usr/bin/perl


# 25/06/2019
# NTM


# Take as arg a filename containing a list of GVCF filenames (with full path, 
# possible gzipped, one file per line) with one or more data columns (ie samples).
# These GVCFs should be produced by Strelka or by this program, if you feed random
# GVCF files you will probably need to adapt the code (although I try to be defensive).
# Produce to stdout a GVCF file, where:
# - header is copied from first file, except: 
#   all INFO descriptions except END and BLOCKAVG_min30p3a are stripped;
#   if --cleanheaders, ##contig headers are stripped except chr1-22 and X,Y,M;
#   ##mergeGVCFs= lines from all infiles are copied;
#   #CHROM line is modified by appending the identifiers of all samples.
# - any line whose FILTER column contains a key from %filtersApplied gets skipped.
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
#   QUAL replaced by '.'
#   INFO and FILTER replaced by '.' except for non-variant blocks:
#      * if all infiles have overlapping blocks with some common FILTERs, we produce a new 
#        block covering the intersection (POS and END get updated, BLOCKAVG_* stays) with
#        the largest common set of FILTER values;
#      * otherwise we print individual lines, not a block, and each sample gets its FT
#        with all its FILTERs;
#   FORMAT contains the union of FORMAT keys from all files (see below);
#   DATA columns contain the CORRECTED sample data for every sample.
#   CORRECTED means:
#     if a key was missing for a sample it gets '.'
#     GT gets corrected values for ALTs
#     GQ, GQX, DPI, DP, DPF, SB, FT, PS, PGT, PID don't change
#     AD/ADF/ADR, get 0 for new ALTs
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
# or same 2 but replace DPI with DP:DPF and add SB:
# GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL
# GT:GQ:GQX:DP:DPF:AD:ADF:ADR:SB:FT:PL:PS
#
# From GATK 4.1.8.1, in VCFs from GenotypeGVCFs via GenomicsDB we have only:
# GT:AD:DP:GQ:PGT:PID:PL:PS
# GT:AD:DP:GQ:PL
# in GATK GVCFs straight from HaplotypeCaller we have:
# GT
# GT:AD:DP:GQ:PGT:PID:PL:PS:SB
# GT:AD:DP:GQ:PL:SB
# GT:DP:GQ:MIN_DP:PL
# GT:GQ:PL
# GT:PGT:PID:PS
#
# In my output files I can change the order of fields and
# create different combinations. But this script can read 
# it's own output.
###########
# NOTES on FILTERS (Strelka2 v2.9.10):
## frequent and seem all relevant: LowGQX LowDepth HighDPFRatio
## never seen without LowGQX or LowDepth: SiteConflict
##  rare, never seen with END, not the most relevant IMO: HighSNVSB
## never seen in my data: IndelConflict NotGenotyped PloidyConflict
#
###########
# NOTES on normalization and merging of variants:
# - as stated, I do NOT left-align variants, it's too much overhead. 
#   I am hoping that Strelka does it correctly, or at least I'm hoping
#   it is consistent in its behavior. If it is consistent, merging
#   will work fine and variants can be left-aligned later if needed.
# - I also do NOT merge lines with different POS values, even if
#   the REFs overlap. I'm afraid doing this would result in large
#   REFs and huge lists of complex ALTs (most of which just describe
#   a SNV...).



use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use Parallel::ForkManager;


#############################################
## hard-coded stuff that shouldn't change much

# max number of lines to read in a single batch (from the first infile,
# number of lines from other infiles will be similar). Each batch is then 
# processed by a worker thread.
# Increasing $batchSize increases the RAM consumption linearly, but should
# improve performance due to better read performance (because better OS buffering),
# less overhead from sub calls, process creations, tmpFiles creations, etc...
# With ~450 samples, each thread uses ~1GB RAM with batchSize==10k in my hands;
# while with 10 samples we're down to 32MB RAM per thread (still batchSize==10k).
my $batchSize = 100000;


# filters to apply: any line whose FILTER value contains a key of %filtersApplied
# will be skipped.
# For performance reasons we do NOT check the FT fields for these keys, we assume
# that the infiles are single-sample Strelka GVCFs or were produced by this script
# with the same %filtersApplied.
# [see NOTES on FILTERS above for how I chose these]
my %filtersApplied = ('LowGQX'=>1, 'LowDepth'=>1, 'HighDPFRatio'=>1);


#############################################
## options / params from the command-line

my $USAGE = '
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--filelist string [required] : file containing a list of GVCF filenames to merge,
    including paths, one file per line
--tmpdir string [default = tmpdir_mergeGVCFs/] : subdir where tmp files will
    be created, must not pre-exist and will be removed after execution
--jobs N [default = 16] : number of parallel jobs=threads to run
--cleanheaders : don\'t print ##contig headers except for chr1-22 and X,Y,M
--help : print this USAGE';


# file containing list of GVCFs, no default
my $fileList;

# for multi-threading, need to create a tmpDir. It will
# be removed when we are done and must not pre-exist.
# to improve performance if you have ample RAM you could
# create a RAMDISK and provide a subdir in it
my $tmpDir = "tmpdir_mergeGVCFs/";

# number of jobs
my $numJobs = 16;

# if $cleanHeaders, don't print ##contig headers except for regular 
# chroms (1-22, X, Y, M)
my $cleanHeaders = '';

# help: if true just print $USAGE and exit
my $help = '';

GetOptions ("filelist=s" => \$fileList,
	    "tmpdir=s" => \$tmpDir,
	    "jobs=i" => \$numJobs,
	    "cleanheaders" => \$cleanHeaders,
	    "help" => \$help)
    or die("E: Error in command line arguments\n\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) &&
    die "$USAGE\n\n";

(-e $tmpDir) && 
    die "E: tmpDir $tmpDir exists, please remove or rename it, or provide a different one with --tmpdir\n";

($fileList) ||
    die "E: you MUST provide a list of GVCFs to merge in a file, with --filelist.\n$USAGE\n";
(-f $fileList) ||
    die "E: provided --filelist $fileList doesn't exist or can't be read\n";


#############################################
## fill @infiles

# array of filehandles, open for reading each GVCF infile
my @infiles = ();

open(FILES, $fileList) ||
    die "E: cannot open provided argument $fileList\n";

while(my $file = <FILES>) {
    chomp($file);
    # skip blank lines
    ($file =~ /^\s*$/) && next;
    (-f $file) || die "E: infile $file from fileList $fileList doesn't exist\n";
    if ($file =~ /\.vcf$/) {
	open(my $infile, $file) || die "E: cannot open infile $file for reading\n";
	push(@infiles, $infile);
    }
    elsif ($file =~ /\.vcf\.gz$/) {
	open(my $infile, "gunzip -c $file |") || die "E: cannot gunzip-open infile $file for reading\n";
	push(@infiles, $infile);
    }
    else {
	die "E: infile $file from fileList $fileList does not end in .vcf or .vcf.gz, why?\n";
    }
}

close(FILES);

# mkdir only now so any previous error doesn't leave an existing tmpDir 
mkdir($tmpDir) || 
    die "E: cannot mkdir tmpDir $tmpDir\n";


#############################################
# deal with headers and fill @numSamples

my $now = strftime("%F %T", localtime);
warn("I: $now - $0 STARTING TO WORK\n");

# same indexes as @infiles, value is the number of samples in the infile
my @numSamples;

# end of header == current command-line and #CHROM line . This allows to
# print the ##mergeGVCFs lines from all files.
my $headerEnd = "";

# copy header from first file
my $infile = $infiles[0];
while(my $line = <$infile>) {
    if (($line =~ /^##INFO/) && ($line !~ /^##INFO=<ID=END,/) && ($line !~ /^##INFO=<ID=BLOCKAVG_min30p3a,/)) {
	next;
    }
    elsif (($cleanHeaders) && ($line =~ /^##contig=<ID=([^,]+),/)) {
	my $contig = $1;
	($contig =~ /^chr[\dXYM]\d?$/) && (print $line);
	# otherwise this is not a regular chrom, don't print
    }
    elsif ($line =~ /^##/) {
	print $line;
    }
    elsif ($line =~ /^#CHROM/) {
	# add full command-line to headers
	my $com = qx/ps -o args= $$/;
	chomp($com);
	#$com .= " < ".`readlink -f /proc/$$/fd/0` ;
	#chomp($com);
	$com .= " > ".`readlink -f /proc/$$/fd/1` ;
	chomp($com);
	$com .= " 2> ".`readlink -f /proc/$$/fd/2` ;
	chomp($com);
	$headerEnd = "##mergeGVCFs=<commandLine=\"$com\">\n";
	chomp($line);
	$headerEnd .= $line;
	# VCF has 9 columns in addition to the data columns
	my @f = split(/\t/,$line);
	push(@numSamples, scalar(@f) - 9);
	last;
    }
    else {
	die "E: parsing header from first infile of $ARGV[0], found bad line:\n$line";
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
	    # grab sample id
	    chomp($line);
	    my @f = split(/\t/,$line);
	    push(@numSamples, scalar(@f) - 9);
	    ($line =~ /\tFORMAT(\t.+)$/) ||
		die "E parsing CHROM header line from file $i:\n$line\n";
	    $headerEnd .= $1;
	    last;
	}
	else {
	    die "E: parsing header from infile $ARGV[$i], found bad line:\n$line";
	}
    }
}

print "$headerEnd\n";

# sanity
(@numSamples == @infiles) || 
    die "numSamples and infiles disagree...\n@numSamples\n@infiles\n";

# flush stdout before starting our eatTmpFiles job
STDOUT->flush();

#############################################
# done with headers, now deal with bodies...

# create fork manager
my $pm = new Parallel::ForkManager($numJobs);

# spawn a child process that waits for workers to finish producing batches,
# and prints the tmpfiles to stdout in correct order, cleaning up behind itself
if (! $pm->start) {
    &eatTmpFiles($tmpDir);
    $pm->finish;
}

# array of "lines", a line is an arrayref to a tab-splitted line from infile, one
# line (max) per infile, each line has been parsed from the infile
#  but didn't belong to the current batch and will be dealt with 
# in the next batch
my @startNextBatch;

# boolean flag, true iff current batch is the last
my $lastBatch = 0;
# initialize with first non-filtered dataline from each file
foreach  my $i (0..$#infiles) {
    my $infile = $infiles[$i];
    my $lineR = &grabNextLine($infile);
    # every infile must have at least one non-filtered dataline
    ($lineR) || 
	die "E: infile $i doesn't have ANY non-filtered dataline?? Something must be wrong";
    push(@startNextBatch, $lineR);
}

# need $batchNum for the tmpFiles
my $batchNum = 0;

while (!$lastBatch) {
    # precondition: if (!$lastBatch) we must have a line
    # for file 0 in @startNextBatch
    $batchNum++;
    # chrom to deal with in the current batch, grab it from first file
    # (we assume first file has at least one non-filtered line with each
    # chrom, this is checked at the end)
    my $thisChr = $startNextBatch[0]->[0] ;

    # array of refs (one per infile) to arrays of (lines == arrayrefs)
    my @batchToMerge;
    # move stored line from first file to @batchToMerge
    $batchToMerge[0] = [$startNextBatch[0]];
    $startNextBatch[0] = undef;

    # smallest position to deal with in next batch, or 0 if all remaining
    # lines for current chrom must be parsed (ie next batch is another chrom or 
    # this is the last batch). Initialize with some non-zero value.
    my $posNextBatch = 1;

    # start with first file, so we can set $posNextBatch
    my $infile = $infiles[0];
    foreach my $j (1..$batchSize) {
	if (my $lineR = &grabNextLine($infile)) {
	    my $chr = $lineR->[0];
	    if ($chr eq $thisChr) {
		&addLineToBatch($batchToMerge[0], $lineR);
	    }
	    else {
		# $lineR is from another chrom, will be for next batch
		$posNextBatch = 0;
		$startNextBatch[0] = $lineR;
		last;
	    }
	}
	else {
	    # no more non-filtered lines in $infile
	    $lastBatch = 1;
	    $posNextBatch = 0;
	    last;
	}
    }
    # 3 possibilities:
    # - we ran out of lines -> $posNextBatch==0
    # - we read a line with a different chrom -> $posNextBatch==0 again
    # - otherwise we have one too many lines in @{$batchToMerge[0]}, move it
    #   to @startNextBatch and set $posNextBatch
    if ($posNextBatch) {
	my $lineR = pop(@{$batchToMerge[0]});
	$posNextBatch = $lineR->[1];
	$startNextBatch[0] = $lineR;
    }

    # Now deal with all files except the first:
    # move relevant previously stored lines of non-first files to @batchToMerge
    foreach my $i (1..$#infiles) {
	if (my $lineR = $startNextBatch[$i]) {
	    my $chr = $lineR->[0];
	    my $pos = $lineR->[1];
	    if (($chr eq $thisChr) && (($pos < $posNextBatch) || ($posNextBatch==0))) {
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
    foreach my $i (1..$#infiles) {
	# if we already have a line in startNextBatch, this file doesn't
	# have any more lines for $thisChr before $posNextBatch
	($startNextBatch[$i]) && next;
	# otherwise:
	my $infile = $infiles[$i];
	while (my $lineR = &grabNextLine($infile)) {
	    my $chr = $lineR->[0];
	    my $pos = $lineR->[1];
	    if (($chr eq $thisChr) && (($pos < $posNextBatch) || ($posNextBatch==0))) {
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

    # OK we can merge the batch and print the result to a tmpfile
    $pm->start && next;
    # if you change $tmpOut or $tmpOutFlag you MUST EDIT &eatTmpFiles()
    my $tmpOut = "$tmpDir/$batchNum.g.vcf";
    open(my $outFH, "> $tmpOut") || die "cannot open $tmpOut for writing\n";
    &mergeBatchOfLines(\@batchToMerge, \@numSamples, $outFH);
    close($outFH);
    # done building this tmpfile, create flag-file
    my $tmpOutFlag = $tmpOut."_done";
    open(OUTFLAG, "> $tmpOutFlag") || die "cannot open flagfile $tmpOutFlag for writing\n";
    print OUTFLAG "$batchNum\n";
    close(OUTFLAG);
    $pm->finish;
}

# some children are still processing batches, but we know the last
# batchNum that will ever exist, tell &eatTmpFiles() so it can exit
# (of course if you change $tmpOutLast you have to edit &eatTmpFiles)
my $tmpOutLast = "$tmpDir/lastBatch";
open(OUTLAST, "> $tmpOutLast") || die "cannot open tmp-last-file $tmpOutLast for writing\n";
print OUTLAST "$batchNum\n";
close OUTLAST;

$pm->wait_all_children;
# sanity: 
foreach my $i (0..$#infiles) {
    ($startNextBatch[$i]) && 
	die "All done but startNextBatch not empty for file $i:\n".join("\t",@{$startNextBatch[$i]})."\n";
}
foreach my $infile (@infiles) {
    close($infile);
}
# clean up: remove $tmpDir 
$now = strftime("%F %T", localtime);

rmdir($tmpDir) || 
    die "E: $now - all done but cannot rmdir tmpDir $tmpDir, why? $!\n";
warn("I: $now - $0 ALL DONE, COMPLETED SUCCESSFULLY\n");



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
# DATA columns, with INFO cleared ('.') except if line is a non-var block
# (ie has END=).
# The last array element has all DATA columns in a single string.
sub grabNextLine {
    (@_ == 1) || die "grabNextLine needs 1 arg.\n";
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
	# if we get here no filters apply, clear INFO if not in a non-var block
	($line[7] =~ /^END=/) || ($line[7] = ".");

	# normalize variants:
	# 1. if length >= 2 for REF and all ALTS, and if REF and all ALTs have 
	#    common ending bases, remove them (keeping at least 1 base everywhere).
	while ($line[3] =~ /\w(\w)$/) {
	    # ref has at least 2 chars
	    my $lastRef = $1;
	    # split here because most lines won't even enter the while loop
	    my @alts = split(/,/,$line[4]);
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
		($line[3] =~ s/$lastRef$//) || 
		    die "WTF can't remove $lastRef from end of $line[3]\n";
		foreach my $i (0..$#alts) {
		    ($alts[$i] =~ s/$lastRef$//) || 
			die "WTF can't remove $lastRef from end of alt $i == $alts[$i]\n";
		}
		$line[4] = join(',',@alts);
	    }
	    else {
		# can't remove $lastRef, get out of while loop
		last;
	    }
	}

	# 2. if length >= 2 for REF and all ALTS, and if REF and all ALTs have 
	#    common starting bases, remove them (keeping at least 1 base everywhere)
	#    and adjust POS.
	while ($line[3] =~ /^(\w)\w/) {
	    my $firstRef = $1;
	    my @alts = split(/,/,$line[4]);
	    my $removeFirst = 1;
	    foreach my $alt (@alts) {
		if ($alt !~ /^$firstRef\w/) {
		    $removeFirst = 0;
		    last;
		}
	    }
	    if ($removeFirst) {
		($line[3] =~ s/^$firstRef//) || 
		    die "WTF can't remove $firstRef from start of $line[3]\n";
		foreach my $i (0..$#alts) {
		    ($alts[$i] =~ s/^$firstRef//) || 
			die "WTF can't remove $firstRef from start of alt $i == $alts[$i]\n";
		}
		$line[4] = join(',',@alts);
		$line[1]++;
	    }
	    else { last; }
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
    (@_ == 2) || die "addLineToBatch needs 2 args.\n";
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
	if ($lineR->[4] eq '.') {
	    # if lineR is non-var ignore it at POS, but make it start at POS+1
	    # if it's a block (just skip it if it's not a block)
	    $lineR = &splitBlockLine($lineR,$lineR->[1]);
	    ($lineR) && (push(@$batchR,$lineR));
	}
	elsif($prevLineR->[4] eq '.') {
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
	    #warn "W: addLineToBatch, 2 lines with ALTs and same POS, ignoring the second line:\n".
	    # join("\t",@$prevLineR)."\n".join("\t",@$lineR)."\n\n";
	}
    }
    else {
	# newPOS < prevPOS, we can switch the 2 lines and recurse,
	# but how far back do we have to go?
	# easy to code with recursive function, but I sure hope
	# this doesn't happen too often!
	warn "I: addLineToBatch is recursing, around $lineR->[0] $lineR->[1]\n";
	pop(@$batchR);
	&addLineToBatch($batchR, $lineR);
	&addLineToBatch($batchR, $prevLineR);
    }
}


###############
# splitBlockLine: args are:
# - $lineR, a ref to an arrayline;
# - $newEnd, an int which must be >= $lineR->POS.
# -> if $lineR doesn't contain a non-variant block: NOOP & return undef
# -> elsif $newEnd < 0: NOOP & return undef
# -> elsif $newEnd >= $lineR->END: NOOP & return undef
# -> else return a ref to a freshly created arrayline, with same 
#    content as $lineR except POS becomes $newEnd+1 and REF becomes 'N';
#    and also modify $lineR: END becomes $newEnd.
sub splitBlockLine {
    (@_ == 2) || die "splitBlockLine needs 2 args.\n";
    my ($lineR,$newEnd) = @_;

    # if $newEnd < 0, NOOP
    ($newEnd >= 0) || (return());
    # if $lineR isn't a block, NOOP
    ($lineR->[7] =~ /^END=(\d+);/) || (return());
    my $oldEnd = $1;
    # if $lineR is a block but ends before or at $newEnd, NOOP
    ($newEnd < $oldEnd) ||  (return());

    # sanity: can't have newEnd positive but < POS
    ($newEnd < $lineR->[1]) &&
	die "E splitBlockLine, newEnd $newEnd smaller than POS in:\n".join("\t",@$lineR)."\n";

    # else make fresh line, initially copy content of $lineR
    my @newLine = @$lineR ;
    # replace POS with $newEnd+1 and REF with N
    $newLine[1] = $newEnd + 1;
    $newLine[3] = 'N';
    # modify END in $lineR
    ($lineR->[7] =~ s/^END=$oldEnd;/END=$newEnd;/) ||
	die "E splitBlockLine, cannot subst END $oldEnd with $newEnd in:\n".join("\t",@$lineR)."\n";
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
    (@_ == 3) || die "mergeBatchOfLines needs 3 args.\n";
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
		if ($fieldsR->[7] =~ /^END=(\d+);/) {
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
		if ($fieldsR->[7] =~ /^END=(\d+);/) {
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
		die "E: in mergeBatchOfLines I want to merge non-var blocks at pos $smallestPos but I have several formats:".
		keys(%seenFormats);
	    print($outFH &mergeLinesNonVarBlock(\@nextToMerge,$numSamplesR, $longestRef));
	}
	else {
	    # build array from %longestFormat, respecting the order specified 
	    # in @maxFormatSorted, and ignoring MIN_DP since we are not in a non-var block
	    my @maxFormatSorted = ('GT','FT','GQ','GQX','DP','DPF','DPI','AD','ADF','ADR','SB','PL','PS','PGT','PID');
	    my @longestFormat = ();
	    foreach my $fkey (@maxFormatSorted) {
		# we must add FT even if it wasn't there
		if (($longestFormat{$fkey}) || ($fkey eq 'FT')) {
		    push(@longestFormat,$fkey);
		}
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
# ref to @longestFormat: array of all FORMAT keys present in at least one line (ignoring MIN_DP)
# Return the string to print, result of merging the lines.
# NOTE: arrays referenced in $toMergeR WILL BE modified (ALTs are changed)
# Preconditions: 
# - MIN_DP must be AFTER DP in FORMAT if it's there (this is checked)
# - any non-variant blockline should have END==POS (but we don't actually
#   examine or use END, so this doesn't matter)
sub mergeLines {
    (@_ == 4) || die "mergeLines needs 4 args.\n";
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
	# if we're in a non-variant position in this file, nothing to do here
	($alts eq '.') && next;
	if ($longestRef ne $ref) {
	    ($longestRef =~ /^$ref(\w+)$/) || 
		die "longestRef $longestRef doesn't start with ref $ref (file $fileIndex), impossible\n".
		join("\t",@{$toMergeR->[$fileIndex]})."\n";
	    my $extraBases = $1;
	    my @fixedAlts = ();
	    foreach my $thisAlt (split(/,/,$alts)) {
		my $fixedAlt = "$thisAlt$extraBases";
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
    # deterministic)
    my @newAlts = sort {(length($a) <=> length($b)) || ($a cmp $b)} keys(%newAlts);
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
	    die "E in mergeLines: longestFormat contains MIN_DP, it MUST NOT!\n";
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
		die "E: in mergeLines MIN_DP is in a FORMAT but before DP! in file $fileIndex:\n".
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
		    #     GQ GQX DPF DPI SB PS PGT PID don't change
		    #     DP gets MIN_DP value if it exists, otherwise doesn't change
		    #     AD ADF ADR get 0 for new ALTs
		    #     FT gets value from FILTER if it didn't exist, otherwise doesn't change
		    #     PL gets correct new values (see explanation in the code)
		    if (!defined $data[$fi]) {
			# data columns can be shorter than FORMAT when last values are unknown
			last;
		    }
		    elsif ($format[$fi] eq "GT") {
			if (($data[$fi] eq '.') || ($data[$fi] eq './.')) {
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
			    die "E: trying to adjust GT $data[$fi] but can't parse it, in file $fileIndex, line:\n".
				join("\t",@{$toMergeR->[$fileIndex]})."\n";
			}
		    }

		    elsif (($format[$fi] eq "GQ") || ($format[$fi] eq "GQX") || 
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

		    elsif (($format[$fi] eq "AD") || ($format[$fi] eq "ADF")|| ($format[$fi] eq "ADR")) {
			my @ADs = split(/,/, $data[$fi]);
			# copy REF AD and remove from @ADs
			$fixedData[$formatIndex{$format[$fi]}] = shift(@ADs);
			# append ADs for all ALTs, inserting '0' except for 
			# alts also present in currentLine
			foreach my $newAltIndex (0..$#newAlts) {
			    if (defined $altsNew2Old[$newAltIndex]) {
				$fixedData[$formatIndex{$format[$fi]}] .= ",".$ADs[$altsNew2Old[$newAltIndex]]; 
			    }
			    else {
				# this newAlt is not in current line, use zero
				$fixedData[$formatIndex{$format[$fi]}] .= ",0";
			    }
			}
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
			die "unknown format key $format[$fi] found, add some code to deal with it!\n".
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

		# if FT is needed but wasn't found, grab data from FILTER
		if ((defined $formatIndex{"FT"}) && ($#fixedData >= $formatIndex{"FT"}) && (! $fixedData[$formatIndex{"FT"}])) {
		    $fixedData[$formatIndex{"FT"}] = $toMergeR->[$fileIndex]->[6];
		}

		# if AD* is needed but wasn't found, use correct number of dummy '.'
		foreach my $fkey ("AD","ADF","ADR") {
		    if ((defined $formatIndex{$fkey}) && ($#fixedData >= $formatIndex{$fkey}) && (! $fixedData[$formatIndex{$fkey}])) {
			# need one value per allele (one REF and all newAlts), use dummy '.'
			$fixedData[$formatIndex{$fkey}] = '.';
			$fixedData[$formatIndex{$fkey}] .= ",." x scalar(@newAlts);
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
# and they all have the same END (this is checked).
# -> if all lines share at least one FILTER value, we return a new blockline 
# using as FILTER the largest set of common FILTER values
# -> else no FILTER value is shared, we return a bunch of single-position lines
# Precondition: FORMAT ends with DP:DPF:MIN_DP (strelka) or DP:GQ:MIN_DP:PL (gatk),
#    this is checked
sub mergeLinesNonVarBlock {
    (@_ == 3) || die "mergeLinesNonVarBlock needs 3 args.\n";
    my ($toMergeR,$numSamplesR, $longestRef) = @_;

    # find first non-null line
    my $firstNonNull = 0;
    while(! $toMergeR->[$firstNonNull]) {
	$firstNonNull++;
    }
    # we can make a new non-var block iff all blocks share some FILTER value(s)
    my $makeNewBlock = 1;
    # shared filters, initialize with filters of first non-null line
    # we use both a ';'-separated string and an array for efficiency
    my $sharedFilters = $toMergeR->[$firstNonNull]->[6];
    my @sharedFilters = split(/;/,$sharedFilters);
    # then update to keep smallest set of shared filters from all files
    foreach my $fileIndex ($firstNonNull+1..$#$toMergeR) {
	($toMergeR->[$fileIndex]) || next;
	my $thisFilter = $toMergeR->[$fileIndex]->[6];
	# quick test, if strings are equal it's fine
	($thisFilter eq $sharedFilters) && next;
	# trick: start at the end so we can splice out any unshared filters
	foreach my $i (reverse(0..$#sharedFilters)) {
	    my $thisShared = $sharedFilters[$i];
	    if ($thisFilter !~ /$thisShared/) {
		# remove thisShared from @sharedFilters
		splice(@sharedFilters,$i,1);
	    }
	}
	# if nothing is shared we can stop now, otherwise update string version
	if (! @sharedFilters) {
	    # incompatible filter values, can't make a new block
	    $makeNewBlock = 0;
	    last;
	}
	else {
	    $sharedFilters = join(';',@sharedFilters);
	}
    }

    # $toPrint will be returned, it can have one or several lines but is always \n-terminated
    my $toPrint = "";
    if ($makeNewBlock) {
	# ok make a new non-variant block: grab start of line from first non-null file,
	# except for REF: use $longestRef because if any infile has non-N it will be there 
	$toPrint = join("\t", @{$toMergeR->[$firstNonNull]}[0..2]);
	$toPrint .= "\t$longestRef";
	$toPrint .= "\t".join("\t", @{$toMergeR->[$firstNonNull]}[4,5]);

	# use smallest set of shared filters
	$toPrint .= "\t$sharedFilters";

	# grab INFO from first non-null file, all files should have the same (this is checked)
	my $info = ${$toMergeR->[$firstNonNull]}[7];
	$toPrint .= "\t$info";

	# grab FORMAT from first non-null file, all files should have the same (this is checked)
	my $format = $toMergeR->[$firstNonNull]->[8];
	$toPrint .= "\t$format";

	# add data columns (even for files that don't have a line at this pos)
	foreach my $fileIndex (0..$#$numSamplesR) {
	    if ($toMergeR->[$fileIndex]) {
		# check INFO
		($toMergeR->[$fileIndex]->[7] eq $info) ||
		    die "E: in mergeLinesNonVarBlock, info $info in first file $firstNonNull is different from ".
		    $toMergeR->[$fileIndex]->[7]." in file $fileIndex\n";
		# check FORMAT
		($toMergeR->[$fileIndex]->[8] eq $format) || 
		    die "E: in mergeLinesNonVarBlock, format $format in first file $firstNonNull is different from ".
		    $toMergeR->[$fileIndex]->[8]." in file $fileIndex\n";

		# print data columns
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
    }

    else {
	# no shared FILTER values, we must put FILTER in data->FT and produce
	# one line per position (using N for REF except for the first line, 
	# otherwise we would need the reference genome in fasta)

	# grab the FILTER value from each file (leave undef if no line in file)
	my @filters = ();
	foreach my $fileIndex (0..$#$toMergeR) {
	    if ($toMergeR->[$fileIndex]) {
		$filters[$fileIndex] = $toMergeR->[$fileIndex]->[6];
	    }
	}

	# use chrom from first non-null file
	my $chrom = $toMergeR->[$firstNonNull]->[0];

	# grab END from first INFO, we will check that all files have the same
	my $infoFirstNonNull = $toMergeR->[$firstNonNull]->[7];
	($infoFirstNonNull =~ /^END=(\d+);/) ||
	    die "E in mergeLinesNonVarBlock: cant grab END from firstNonNull $firstNonNull\n";
	my $end = $1;

	# FORMAT can also be grabbed from the first file, we will check that all
	# files have the same, but we strip MIN_DP and add FT
	my $formatFirstNonNull = $toMergeR->[$firstNonNull]->[8];
	my $format = $formatFirstNonNull;
	# remove MIN_DP
	# behavior depends on $caller: "strelka" or "gatk"
	# if this code changes you MUST also change the substitutions of $data below 
	my $caller;
	if ($format =~ s/:DP:DPF:MIN_DP$/:DP:DPF/) {
	    $caller = "strelka";
	}
	elsif ($format =~ s/:DP:GQ:MIN_DP:PL$/:DP:GQ:PL/) {
	    $caller = "gatk";
	}
	else {
	    die "E in mergeLinesNonVarBlock: no shared filters but cannot remove MIN_DP from FORMAT in first file $firstNonNull, format is $format\n";
	}
	# add FT right after GT
	($format =~ s/^GT:/GT:FT:/) || 
	    die "E in mergeLinesNonVarBlock: no shared filters but cannot add FT to FORMAT in first file $firstNonNull, format is $format\n";

	my $firstPos = $toMergeR->[$firstNonNull]->[1];
	foreach my $pos ($firstPos..$end) {
	    # use N for REF except for first line
	    my $line = "$chrom\t$pos\t.";
	    if ($pos == $firstPos) {
		$line .= "\t$longestRef\t.\t.";
	    }
	    else {
		$line .= "\tN\t.\t.";
	    }
	    # FILTER becomes '.' (will be added to DATA columns)
	    $line .= "\t.";
	    # INFO becomes '.'
	    $line .= "\t.";
	    # FORMAT
	    $line .= "\t$format";
	    # geno-data: remove MIN_DP, replace DP value with MIN_DP value,
	    # and add the correct FT for each, but go all the way to $#$numSamplesR
	    foreach my $fileIndex (0..$#$numSamplesR) {
		if ($toMergeR->[$fileIndex]) {
		    # sanity: check INFO and FORMAT
		    ($toMergeR->[$fileIndex]->[7] eq $infoFirstNonNull) ||
			die "E in mergeLinesNonVarBlock: no shared, pos $pos, INFO differs in file $fileIndex from $infoFirstNonNull\n"; 
		    ($toMergeR->[$fileIndex]->[8] eq $formatFirstNonNull) ||
			die "E in mergeLinesNonVarBlock: no shared, pos $pos, FORMAT differs in file $fileIndex from $formatFirstNonNull\n";
		    my @dataCols = split("\t",$toMergeR->[$fileIndex]->[9]);
		    foreach my $data (@dataCols) {
			# if no data, leave as '.'
			if ($data ne '.') {
			    # remove MIN_DP and replace DP with its value
			    if ($caller eq "strelka") {
				($data =~ s/:\d+:(\d+):(\d+)$/:$2:$1/) ||
				    die "E in mergeLinesNonVarBlock: no shared, caller $caller, cannot DP->MIN_DP in data from file $fileIndex:$data\n";
			    }
			    elsif ($caller eq "gatk") {
				($data =~ s/:\d+:(\d+):(\d+):([^:]+)$/:$2:$1:$3/) ||
				    die "E in mergeLinesNonVarBlock: no shared, caller $caller, cannot DP->MIN_DP in data from file $fileIndex:$data\n";
			    }
			    else {
				die "E in mergeLinesNonVarBlock: caller $caller not implemented!\n";
			    }
			    # add filters in second field, after GT
			    ($data =~ s/^([^:]+):/$1:$filters[$fileIndex]:/) ||
				die "E in mergeLinesNonVarBlock: no shared, cannot add FT values in data from file $fileIndex:$data\n"; 
			}
			$line .= "\t$data";
		    }
		}
		else {
		    # this file doesn't have data here, use "."
		    foreach my $j (1..$numSamplesR->[$fileIndex]) {
			$line .= "\t.";
		    }
		}
	    }
	    $toPrint .= "$line\n";
	}
    }

    return($toPrint);
}




###############
# this function waits for flagfiles to be created and "eats" the 
# corresponding tmpFiles in order, starting at 1.
# "eating" means print to stdout and remove.
# we also watch for $tmpOutLast, a file that will tell us
# the last batch number to wait for
sub eatTmpFiles {
    (@_ == 1) || die "eatTmpFiles needs 1 arg.\n";
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

	if (-e $tmpOutFlag) {
	    # next batch tmp file is ready
	    open (IN, $tmpOut) || 
		die "E in eatTmpFiles, flagfile $tmpOutFlag exists but cant open tmpout $tmpOut: $!\n";
	    while(<IN>) {
		print;
	    }
	    close(IN);
	    (unlink($tmpOut,$tmpOutFlag) == 2) ||
		die "E in eatTmpFiles, done with files for batch $nextBatch but cannot unlink (both of) them: $!\n";

	    my $now = strftime("%F %T", localtime);
	    warn("I: $now - done printing results from batch $nextBatch\n");
	    $nextBatch++;
	    next;
	}

	elsif (-e $tmpOutLast) {
	    open (IN, $tmpOutLast) || 
		die "E in eatTmpFiles, cannot open tmpOutLast $tmpOutLast although it exists: $!\n";
	    $lastBatch = <IN>;
	    chomp($lastBatch);
	    close(IN);
	    unlink($tmpOutLast) || 
		die "E in eatTmpFiles, cannot unlink tmpOutLast $tmpOutLast: $!\n";
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

