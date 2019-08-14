#!/usr/bin/perl


# 05/04/2018
# NTM

# Take as arg a list of GVCF filenames (with full path, possible gzipped) 
# with one or more data columns (ie samples);
# produce to stdout a GVCF file, where:
# - header is copied from first file, except: 
#   all INFO descriptions are stripped
#   SB description in FORMAT is stripped
#   #CHROM line is modified by appending the identifiers of all samples.
# - any chrom:coord line present in any infile will appear in the output, with:
#   ID set to '.';
#   the longest REF from all infiles;
#   all ALTs from all infiles, adjusted to fit the longest REF (append extra 
#      bases if needed);
#   QUAL, FILTER and INFO replaced by '.';
#   FORMAT contains the longest FORMAT string from all files, stripped of SB;
#   DATA columns contain the CORRECTED sample data for every sample, 
#     or './.' if a sample doesn't have anything called for that chrom:coord.
#     CORRECTED means:
#     DP and GQ don't change
#     GT becomes ./.
#     AD and SAC get NONREF values for new ALTs
#     PL gets correct new values (see explanation in the code)
#     SB gets stripped
#     PGT and PID get '.' if missing.
#
# I think I've reverse-engineered CombineGVCFs correcty.
# So, I'll implement this. Then I can make sure I got it right by just comparing
# the TSVs I produce from going through mergeGvcfs.pl vs straight through GATK.


use strict;
use warnings;
use Parallel::ForkManager;

my $numJobs = 12 ;
my $pm = new Parallel::ForkManager($numJobs);

# max number of lines to read in a single batch from the first infile,
# number of lines from other files may vary around this but must
# satisfy the &mergeBatchOfLines preconditions.
# With batchSize==10k and ~400 samples, each thread uses ~3.5GB RAM in my hands.
# Increasing $batchSize increases the RAM consumption linearly, but also
# reduces the number of tmpFiles created, and probably improves read 
# performance (because better OS buffering) and cuts down
# on the overhead from sub calls, process creations, etc...
my $batchSize = 10000;

# for multi-threading, need to create a tmpDir. It will
# be removed when we are done but must not pre-exist.
# MAKE SURE THIS IS DEFINED AND REASONABLE, WE "rm -r" IT AT THE END!!!
# ... NO actually we don't, I'm chickening out. I'll just tell the user to do it...
my $tmpDir = "tmpdir_mergeGVCFs/";
(-e $tmpDir) && die "tmpDir $tmpDir exists, please remove or rename it\n";

(@ARGV == 0) &&
    die "This script must be called with a list of VCF filenames as arguments\n";

# array of filehandles, open for reading each infile
my @infiles = ();

foreach my $file (@ARGV) {
    (-f $file) || die "E: infile $file doesn't exist\n";
    if ($file =~ /\.vcf$/) {
	open(my $infile, $file) || die "E: cannot open infile $file for reading\n";
	push(@infiles, $infile);
    }
    elsif ($file =~ /\.vcf\.gz$/) {
	open(my $infile, "gunzip -c $file |") || die "E: cannot open infile $file for reading\n";
	push(@infiles, $infile);
    }
    else {
	die "argument $file does not end in .vcf or .vcf.gz, why?\n";
    }
}

mkdir($tmpDir) || die "cannot mkdir tmpDir $tmpDir\n";


#############################################
# deal with headers and fill @numSamples

# same indexes as @infiles, value is the number of samples in the infile
my @numSamples;

# copy header from first file
my $infile = $infiles[0];
while(my $line = <$infile>) {
    if (($line =~ /^##INFO/) || ($line =~ /^##FORMAT=<ID=SB,/)) {
	next;
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
	print "##mergeGVCFs=<commandLine=\"$com\">\n";
	chomp($line);
	print $line;
	# VCF has 9 columns in addition to the data columns
	my @f = split(/\t/,$line);
	push(@numSamples, scalar(@f) - 9);
	last;
    }
    else {
	die "E: parsing header from first infile $ARGV[0], found bad line:\n$line";
    }
}

# skip headers from other infiles but grab sample ids
foreach my $i (1..$#infiles) {
    my $infile = $infiles[$i];
    while(my $line = <$infile>) {
	if ($line =~ /^##/) {
	    #NOOP, skip
	}
	elsif ($line =~ /^#CHROM/) {
	    # grab sample id
	    chomp($line);
	    my @f = split(/\t/,$line);
	    push(@numSamples, scalar(@f) - 9);
	    ($line =~ /\tFORMAT(\t.+)$/) ||
		die "E parsing CHROM header line from file $i:\n$line\n";
	    print $1;
	    last;
	}
	else {
	    die "E: parsing header from infile $ARGV[$i], found bad line:\n$line";
	}
    }
}

print "\n";

# sanity
(@numSamples == @infiles) || 
    die "numSamples and infiles disagree...\n@numSamples\n@infiles\n";


#############################################
# done with headers, now deal with bodies...

# array of "lines", a line is a ref to the chomped string from infile, one
# line (max) per infile, each line has been parsed from the infile
#  but didn't belong to the current batch and will be dealt with 
# in the next batch
my @startNextBatch;

# boolean flag, true iff current batch is the last
my $lastBatch = 0;
# initialize with first dataline from each file
foreach  my $i (0..$#infiles) {
    my $infile = $infiles[$i];
    my $line = <$infile>;
    ($line) || die "no dataline in infile $i == $ARGV[$i]\n";
    chomp($line);
    push(@startNextBatch, \$line);
}

# need $batchNum for the tmpFiles
my $batchNum = 0;

while (!$lastBatch) {
    # precondition: if (!$lastBatch) we must have a line
    # for file 0 in @startNextBatch
    $batchNum++;
    # chrom to deal with in the current batch, grab it from first file
    # (we assume first file has at least one line with each chrom,
    # this is checked at the end)
    (${$startNextBatch[0]} =~ /^([^\t]+)\t/) || 
	die "cannot grab chrom from first file in batch $batchNum\n".${$startNextBatch[0]}."\n";
    my $thisChr = $1;

    # array of refs (one per infile) to arrays of lines == stringrefs
    my @batchToMerge;
    # move stored line from first file to @batchToMerge
    $batchToMerge[0] = [$startNextBatch[0]];
    $startNextBatch[0] = undef;

    # smallest position to deal with in next batch, or 0 if all remaining
    # lines for current chrom must be parsed (ie next batch is another chrom or 
    # this is the last batch). Initialize with some non-zero value.
    my $nextPos = 1;

    # start with first file, so we can set $nextPos
    my $infile = $infiles[0];
    foreach my $j (1..$batchSize) {
	if (my $line =  <$infile>) {
	    chomp($line);
	    ($line =~ /^([^\t]+)\t/) || die "cannot grab chrom in:\n$line\n";
	    my $chr = $1;
	    if ($chr eq $thisChr) {
		push(@{$batchToMerge[0]}, \$line);
	    }
	    else {
		# $line is from another chrom, will be for next batch
		$nextPos = 0;
		$startNextBatch[0] = \$line;
		last;
	    }
	}
	else {
	    # no more lines in $infile
	    $lastBatch = 1;
	    $nextPos = 0;
	    last;
	}
    }
    # 3 possibilities:
    # - we ran out of lines -> $nextPos==0
    # - we read a line with a different chrom -> $nextPos==0 again
    # - otherwise we have one too many lines in @linesThisFile, move it
    #   to @startNextBatch and set $nextPos
    if ($nextPos) {
	my $lineR = pop(@{$batchToMerge[0]});
	($$lineR =~ /^[^\t]+\t(\d+)\t/) || die "cannot grab pos in:\n$$lineR\n";
	$nextPos = $1;
	$startNextBatch[0] = $lineR;
    }

    # move relevant previously stored lines to @batchToMerge
    foreach my $i (1..$#infiles) {
	if ($startNextBatch[$i]) {
	    (${$startNextBatch[$i]} =~ /^([^\t]+)\t(\d+)\t/) ||
		die "cannot grab chrom and pos in:\n${$startNextBatch[$i]}\n";
	    my ($chr,$pos) = ($1,$2);
	    if (($chr eq $thisChr) && (($pos < $nextPos) || ($nextPos==0))) {
		$batchToMerge[$i] = [$startNextBatch[$i]];
		$startNextBatch[$i] = undef;
	    }
	    # else: file $i doesn't have any more lines for this chrom before $nextPos,
	    # leave its next line in @startNextBatch ie NOOP
	}
    }


    # Now deal with all files except the first
    foreach my $i (1..$#infiles) {
	# if we already have a line in startNextBatch, this file doesn't
	# have any more lines for $thisChr before $nextPos
	($startNextBatch[$i]) && next;
	# otherwise:
	my $infile = $infiles[$i];
	while (my $line = <$infile>) {
	    chomp($line);
	    ($line =~ /^([^\t]+)\t(\d+)\t/) ||
		die "cannot grab chrom and pos in:\n$line\n";
	    my ($chr,$pos) = ($1,$2);
	    if (($chr eq $thisChr) && (($pos < $nextPos) || ($nextPos==0))) {
		push(@{$batchToMerge[$i]},\$line);
	    }
	    else {
		$startNextBatch[$i] = \$line;
		last; # last line for infile $i
	    }
	}
    }

    # OK we can merge the batch and print the result to a tmpfile
    # the tmpFile names should be left-padded with zeroes so we can cat *
    my $batchNumPadded = "0" x (10-length($batchNum));
    $batchNumPadded .= $batchNum;
    my $tmpOut = "$tmpDir/$batchNumPadded.g.vcf";
    $pm->start && next;
    open(my $outFH, "> $tmpOut") || die "cannot open $tmpOut for writing\n";
    &mergeBatchOfLines(\@batchToMerge, \@numSamples, $outFH);
    close($outFH);
    $pm->finish;
}

$pm->wait_all_children;
# sanity: 
foreach my $i (0..$#infiles) {
    ($startNextBatch[$i]) && 
	die "All done but startNextBatch not empty for file $i:\n".${$startNextBatch[$i]}."\n";
}
foreach my $infile (@infiles) {
    close($infile);
}

system("cat $tmpDir/*");
# I HATE DOING THIS!!
#system("rm -r $tmpDir");
# I'm chickening out: just tell the user to do it
warn("All done, you should clean up the tmpDir by running:\nrm -r $tmpDir\n");

# mergeBatchOfLines: args are
# $batchToMergeR == arrayref of refs (one per infile) to arrays of 
#    "lines", each line is actually a stringref to the chomped line
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
#   chrom:pos in the NEXT batch, must be in this batch.
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
	# array-lines are refs to arrays that result from splitting 
	# the lines on \t and clearing INFO
	my @nextToMerge ;
	# current smallest pos
	my $nextPos = 0;
	# longest REF among the toMerge lines
	my $longestRef = "";
	# FORMAT keys that appear in the toMerge lines
	my %longestFormat = ();
	# FORMAT strings that appear in toMerge lines
	my %seenFormats = ();
	# boolean flag: true iff lines to merge don't have any ALTs
	my $empty;

	foreach my $i (0..$#infiles) {
	    # ignore files without remaining lines in this batch
	    ($batchToMergeR->[$i]->[0]) || next;

	    my @fields = split(/\t/, ${$batchToMergeR->[$i]->[0]});
	    # clear INFO
	    $fields[7] = ".";
	    my ($pos,$ref,$alts,$format) = @fields[1,3,4,8];
	    if (($pos < $nextPos) || (! $nextPos)) {
		# this is the new smallest
		@nextToMerge = ();
		$nextToMerge[$i] = \@fields;
		$nextPos = $pos;
		$longestRef = $ref;
		%longestFormat = ();
		foreach my $fkey (split(/:/,$format)) {
		    $longestFormat{$fkey} = 1;
		}
		%seenFormats = ();
		$seenFormats{$format} = 1;
		$empty = 1;
		($alts eq '<NON_REF>') || ($empty = 0);
	    }
	    elsif ($pos == $nextPos) {
		$nextToMerge[$i] = \@fields;
		(length($ref) > length($longestRef)) && ($longestRef = $ref);
		if (! $seenFormats{$format}) {
		    foreach my $fkey (split(/:/,$format)) {
			$longestFormat{$fkey} = 1;
		    }
		    $seenFormats{$format} = 1;
		}
		($alts eq '<NON_REF>') || ($empty = 0);
	    }
	    # else current file has "larger" position, NOOP
	}

	# remove the chosen lines from $batchToMergeR and decrement 
	# $filesNotFinished if we finished a file
	foreach my $i (0..$#nextToMerge) {
	    if ($nextToMerge[$i]) {
		shift(@{$batchToMergeR->[$i]});
		($batchToMergeR->[$i]->[0]) || ($filesNotFinished--);
	    }
	}

	if ($empty) {
	    print($outFH &mergeLinesAllEmpty(\@nextToMerge,$numSamplesR));
	}
	else {
	    # build array from %longestFormat and remove SB, respecting the 
	    # order specified in @maxFormatSorted (which DOESN'T have SB, this
	    # is how we discard SB)
	    my @maxFormatSorted = ('GT','AD','DP','GQ','PGT','PID','PL','SAC');
	    my @longestFormat = ();
	    foreach my $fkey (@maxFormatSorted) {
		($longestFormat{$fkey}) && push(@longestFormat,$fkey);
	    }
	    print($outFH &mergeLines(\@nextToMerge, $numSamplesR, \$longestRef, \@longestFormat));
	}
    }
}


# mergeLines: merge (up to) one line per infile, args are:
# arrayref of refs to array-lines, one array-line per infile, 
#    all array-lines have same chrom:pos or are undef
# ref to @numSamples == same indexes as @infiles, value is the number of samples in the infile
# ref to $longestRef wich contains the longest REF among the strings
# ref to @longestFormat: array of all FORMAT keys present in at least one line (except SB)
# Return the string to print, result of merging the strings
# NOTE: arrays referenced in $toMergeR WILL BE modified
sub mergeLines {
    (@_ == 4) || die "mergeLines needs 4 args.\n";
    my ($toMergeR, $numSamplesR, $longestRefR, $longestFormatR) = @_;

    # %newAltsCount: this is filled in step 1 and used in step 2 
    # key == a (new, fixed) ALT that exists at this chrom:pos, 
    # value == number of samples where it may exist (ie samples
    # appearing in files where this ALT appears).
    # NONREF is not used as key but we know it exists (and must go last).
    my %newAltsCount;

    ####################################
    # STEP 1: fix ALTs in @$toMergeR so they correspond 
    # to a REF == $longestRef, and fill %newAltsCount
    foreach my $fileIndex (0..$#$toMergeR) {
	($toMergeR->[$fileIndex]) || next;
	my ($ref,$alts) = @{$toMergeR->[$fileIndex]}[3,4];
	if ($$longestRefR ne $ref) {
	    ($$longestRefR =~ /^$ref(\w+)$/) || 
		die "longestRef $$longestRefR doesn't start with ref $ref (file $fileIndex), impossible\n".
		join("\t",@{$toMergeR->[$fileIndex]})."\n";
	    my $extraBases = $1;
	    my @fixedAlts;
	    foreach my $thisAlt (split(/,/,$alts)) {
		if ($thisAlt eq '<NON_REF>') {
		    # never change NONREF
		    push(@fixedAlts, $thisAlt);
		}
		else {
		    my $fixedAlt = "$thisAlt$extraBases";
		    push(@fixedAlts, $fixedAlt);
		    ($newAltsCount{$fixedAlt}) || ($newAltsCount{$fixedAlt} = 0);
		    $newAltsCount{$fixedAlt} += $$numSamplesR[$fileIndex];
		}
	    }
	    # now fix the actual line: replace ALT
	    my $fixedAltsString = join(',',@fixedAlts);
	    $toMergeR->[$fileIndex]->[4] = $fixedAltsString;
	}
	else {
	    # else $ref is the longestRef, just update newAltsCount
	    foreach my $thisAlt (split(/,/,$alts)) {
		($thisAlt eq '<NON_REF>') && next;
		($newAltsCount{$thisAlt}) || ($newAltsCount{$thisAlt} = 0);
		$newAltsCount{$thisAlt} += $$numSamplesR[$fileIndex];
	    }
	}
    }

    ####################################
    # STEP 2 (NO LOOP): build @newAlts and %newAlts, sorted appropriately
    my @newAlts = sort {$newAltsCount{$b} <=> $newAltsCount{$a}} keys(%newAltsCount);
    push(@newAlts, '<NON_REF>');
    # also fill %newAlts: key==ALT, value==index of that alt in @newAlts
    my %newAlts;
    foreach my $i (0..$#newAlts) {
	$newAlts{$newAlts[$i]} = $i;
    }

    ####################################
    # STEP 3 (NO LOOP): start building line that will be returned
    my $toPrint;
    # grab CHROM POS from first non-null line
    my $firstNonNull = 0;
    while(! $toMergeR->[$firstNonNull]) {
	$firstNonNull++;
    }
    $toPrint = $toMergeR->[$firstNonNull]->[0]."\t".$toMergeR->[$firstNonNull]->[1];
    # ID REF
    $toPrint .= "\t.\t".$$longestRefR;
    # ALT
    $toPrint .= "\t".join(',',@newAlts);
    # QUAL FILTER INFO FORMAT
    $toPrint .= "\t.\t.\t.\t".join(':',@$longestFormatR);


    ####################################
    # STEP 4: fix the data columns and add them to $toPrint
    foreach my $fileIndex (0..$#$numSamplesR) {
	if (! $toMergeR->[$fileIndex]) {
	    # no line, use blank data for all samples from this file
	    foreach my $j (1..$$numSamplesR[$fileIndex]) {
		$toPrint .= "\t./.";
	    }
	}
	else {
	    my @alts = split(/,/,$toMergeR->[$fileIndex]->[4]);
	    # @altsNew2Old: same indexes as @newAlts,
	    # value is the index of $newAlts[$i] in @alts, or undef if it's not there
	    my @altsNew2Old;
	    foreach my $i (0..$#alts) {
		$altsNew2Old[$newAlts{$alts[$i]}] = $i;
	    }

	    my @format = split(/:/,$toMergeR->[$fileIndex]->[8]);

	    # deal with each DATA column
	    foreach my $j (1..$$numSamplesR[$fileIndex]) {
		my @data = split(/:/, $toMergeR->[$fileIndex]->[8+$j]);
		# we will make sure fields in @format are in the same order as
		# in @$longestFormatR (except SB, which we will skip)

		# fixed DATA string for this sample
		my $fixedData;
		# index in @format of next element to look at
		my $fi = 0;
		foreach my $lfi (0..$#$longestFormatR) {
		    #     GT becomes ./.
		    #     DP GQ PGT and PID don't change
		    #     AD and SAC get NONREF values for new ALTs
		    #     PL gets correct new values (see explanation in the code)
		    #     SB gets stripped
		    if (($fi <= $#format) && ($format[$fi] eq $longestFormatR->[$lfi])) {
			if ($format[$fi] eq "GT") {
			    # VCF spec says GT must be first, so no .= or ':'
			    $fixedData = "./.";
			}

			elsif (($format[$fi] eq "DP") || ($format[$fi] eq "GQ") || 
			       ($format[$fi] eq "PGT") || ($format[$fi] eq "PID")) {
			    # just copy values
			    $fixedData .= ":".$data[$fi];
			}

			elsif ($format[$fi] eq "AD") {
			    my @ADs = split(/,/, $data[$fi]);
			    # copy REF AD and remove from @ADs
			    $fixedData .= ":".shift(@ADs);
			    # append ADs for all ALTs, inserting the NON_REF value except for 
			    # alts also present in currentLine
			    foreach my $newAltIndex (0..$#newAlts) {
				if (defined $altsNew2Old[$newAltIndex]) {
				    $fixedData .= ",".$ADs[$altsNew2Old[$newAltIndex]]; 
				}
				else {
				    # this newAlt is not in current line, use NONREF value
				    $fixedData .= ",".$ADs[$#ADs];
				}
			    }
			}

			elsif ($format[$fi] eq "SAC") {
			    my @SACs = split(/,/, $data[$fi]);
			    # grab REF SACs and remove from @SACs
			    $fixedData .= ":".shift(@SACs);
			    $fixedData .= ",".shift(@SACs);
			    # fill SACs for all ALTs, inserting NON_REF SACs for absent alts
			    foreach my $newAltIndex (0..$#newAlts) {
				if (defined $altsNew2Old[$newAltIndex]) {
				    $fixedData .= ",".$SACs[2 * $altsNew2Old[$newAltIndex]];
				    $fixedData .= ",".$SACs[2 * $altsNew2Old[$newAltIndex] + 1];
				}
				else {
				    # this newAlt is not in current line, use NONREF values
				    $fixedData .= ",".$SACs[$#SACs - 1].",".$SACs[$#SACs];
				}
			    }
			}

			elsif ($format[$fi] eq "PL") {
			    # ALTs that didn't exist got the NON_REF counts, so
			    # the corresponding PLs will also get the NON_REF values
			    # The VCF spec says the PL for x/y genotype is at index x + y*(y+1)/2
			    # Note that x and y here start at 1 for ALTs, 0 is the REF.
			    # When y is an unknown ALT:
			    # the 0/y PL is copied from 0/NONREF;
			    # the x/y or y/x PL (where x is a previously existing ALT) is 
			    # copied from x/NONREF;
			    # and the y/y or Y1/Y2 PL is copied from NONREF/NONREF == $PLs[$#PLs].
			    # SO in fact the new PL for X/Y is:
			    # 1. replace X by NONREF if X didn't exist before, same with Y;
			    # 2. find previous PL of X/Y.
			    my @PLs = split(/,/, $data[$fi]);

			    my @newPLs;

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
					# X is REF -> NOOP
				    }
				    elsif (defined $altsNew2Old[$x]) {
					$oldX = $altsNew2Old[$x] + 1;
				    }
				    else {
					# X unknown, use index of NONREF in @alts
					$oldX = $#alts + 1;
				    }

				    if ($y == -1) {
					# Y is REF -> NOOP
				    }
				    elsif (defined $altsNew2Old[$y]) {
					$oldY = $altsNew2Old[$y] + 1;
				    }
				    else {
					# Y unknown, use index of NONREF in @alts
					$oldY = $#alts + 1;
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
				    $newPLs[$plDest] = $PLs[$plSource];
				}
			    }
			    $fixedData .= ":".join(',', @newPLs);
			}

			else {
			    die "unknown format key $format[$fi] found, add some code to deal with it!\n".
				join("\t",@{$toMergeR->[$fileIndex]})."\n";
			}

			# in any case we ate up a key from @format and @data, increment $fi
			$fi++;
			# if next field is SB, rip it out
			(defined $format[$fi]) && ($format[$fi] eq 'SB') && ($fi++);
		    }

		    else {
			# $lfi format key is not in @format:
			# use dummy value, hoping any missing fields are single-valued
			$fixedData .= ":.";
		    }
		}

		# make sure every field from @format was dealt with
		if ($fi != $#format + 1) {
		    warn "line has format field at index $fi not present in longestFormat ".
			join(':',@$longestFormatR).
			" or in wrong order, any fields from index $fi onwards were ignored:\n".
			join("\t",@{$toMergeR->[$fileIndex]})."\n";
		}
		# AOK, remove trailing empty fields and print fixed data
		while ($fixedData =~ s/:.$//) {
		    # NOOP
		}
		$toPrint .= "\t$fixedData";
	    }
	}
    }

    ####################################
    # STEP 5: all done! return result
    return("$toPrint\n");
}



# special case of mergeLines: merge lines where the only ALT allele 
# in every non-null line is NON_REF
sub mergeLinesAllEmpty {
    (@_ == 2) || die "mergeLinesAllEmpty needs 2 args.\n";
    my ($toMergeR,$numSamplesR) = @_;
    # grab all columns except DATA from first non-null line
    my $firstNonNull = 0;
    while(! $toMergeR->[$firstNonNull]) {
	$firstNonNull++;
    }
    my $toPrint = join("\t", @{$toMergeR->[$firstNonNull]}[0..8]);

    # add data columns
    foreach my $fileIndex (0..$#$numSamplesR) {
	if ($toMergeR->[$fileIndex]) {
	    foreach my $j (1..$$numSamplesR[$fileIndex]) {
		# replace GT with ./.
		my $data = $toMergeR->[$fileIndex]->[8+$j];
		$data =~ s~^0/0:~./.:~;
		$toPrint .= "\t$data";
	    }
	}
	else {
	    # file $fileIndex doesn't have a line for this position, print 
	    # "blank" data for each sample from infile
	    foreach my $j (1..$$numSamplesR[$fileIndex]) {
		$toPrint .= "\t./.";
	    }
	}
    }

    return("$toPrint\n");
}
