#!/usr/bin/perl

# 13/06/2019
# NTM


# Sometimes we get several pairs of FASTQ files for a single sample,
# eg when the sample was sequenced on several lanes.
# Also, the naming conventions are always different.
#
# This script creates a single pair of "files" (symlinks when possible)
# in $inPath/../$outDir/ for each grexome sample, with uniform naming scheme
# grexomeXXXX_[12].fq.gz .
#
# Basically we want to symlink when there's a single pair of files,
# and otherwise cat @infiles1 > $outDir/$grexome_1.fq.gz
# and same for @infiles2 (in same order).
# NOTE: we CAN just concatenate multiple gzip files and get a sane gzip file,
# see gzip spec https://tools.ietf.org/html/rfc1952 (grep for "members")
#
# We try to find original FASTQs automatically, by looking for .gz files
# in subdirs of $inPath whose name contains specimenID from the xlsx file.
# This is a bit fragile but we sanity-check as much as possible: we make
# sure each source file is used at most once, and we report at the end
# if every source fastq was actually used. Always do a dry run initially
# and examine the log! Then add --real to do the real work.
#
# See $USAGE for arguments.

use strict;
use warnings;
use Getopt::Long;
use Spreadsheet::XLSX;

# make glob case-insensitive (for eg */P17if085*)
use File::Glob qw(:globally :nocase);


#############################################
## options / params from the command-line

my $USAGE = '
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--grexome2sample fileWithPath [required] : full path to an xlsx file
    with "sampleID" in some col and "specimenID" in another 
    (eg .../patient_summary.xlsx). 
--inpath path : path to a directory containing subdirs, each containing the
    actual original fastq.gz files (eg "FASTQs_grexome/").
--outdir subdir : name of brother-dir of $inPath where grexome*.fq.gz files are 
    symlinked/produced [default: "FASTQs_All_Grexomized"]. NOTE: outDir cannot 
    contain a slash, it will be created if needed and will be a BROTHER of the 
    final subdir of $inPath ie $inPath/../$outDir.
--real : actually do the work, default without this is a dry run ie it just
    prints info on what would be done.
--help : print this USAGE';


# filename with path to patient_summary*.xlsx file
my $grex2sampleFile = '';

# path containing subdirs containing the original fq.gz files
my $inPath = '';

# brother-dir of $inPath where new FASTQs are created/symlinked, must 
# be a plain dirname (no slashes allowed)
my $outDir = "FASTQs_All_Grexomized";

# real: if true actually process files, otherwise just print INFO 
# messages on what would be done. Default to false (==dry run)
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

GetOptions ("grexome2sample=s" => \$grex2sampleFile,
	    "inpath=s" => \$inPath,
	    "outdir=s" => \$outDir,
	    "real" => \$real,
	    "help" => \$help)
    or die("Error in command line arguments\n\n$USAGE\n");


# make sure required options were provided and sanity check them
($help) &&
    die "$USAGE\n\n";
($grex2sampleFile) || 
    die "$USAGE\n\nE: you MUST specify the path&filename of the patient_summary*.xlsx file, with --grexome2sample\n";
(-f $grex2sampleFile) ||
    die "E: the supplied grexome2sampleFile doesn't exist\n";
($inPath) || 
    die "$USAGE\n\nE: you MUST specify the path holding subdirs that contain the original FASTQs, with --inpath\n";
(-d $inPath) ||
    die "E: the provided inpath $inPath doesn't exist as a dir\n";
($outDir =~ m~/~) &&
    die "$USAGE\n\nE: the provided outdir $outDir has a slash, it must be a plain name\n";

# strip any trailing slashes from $inPath, and split $inPath into $parentDir and $lastDir
$inPath =~ s~/*$~~ ;
my ($parentDir,$lastDir);
if ($inPath =~ m~^(.+/)([^/]+)$~) {
    ($parentDir,$lastDir) = ($1,$2);
}
else {
    ($parentDir,$lastDir) = ("",$inPath);
}

# create $outDir if needed
(-d "$parentDir$outDir/") || (mkdir("$parentDir$outDir/")) || 
    die "E: outDir $parentDir$outDir doesn't exist as a dir but can't be created\n";


#############################################
## parse $grex2sampleFile and fill @grex2sample

# $grex2sample[$grexNum] will hold the "specimen" value from the xlsx file,
# corresponding to grexome number $grexNum
# EXCEPTION for Strasbourg (grabbed from "Center col) where we prepend 'PRY-'
# because it is present in the filenames and the specimen names
# are so stupidly short they are ambiguous
my @grex2sample = ();
{
    my $workbook = Spreadsheet::XLSX->new("$grex2sampleFile");
    (defined $workbook) ||
	die "E when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($grexCol, $specimenCol, $centerCol) = (-1,-1,-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# skip empty columns
	($cell) || next;
	($cell->value() eq "sampleID") &&
	    ($grexCol = $col);
	($cell->value() eq "specimenID") &&
	    ($specimenCol = $col);
	($cell->value() eq "Center") &&
	    ($centerCol = $col);
    }
    ($grexCol >= 0) ||
	die "E parsing xlsx: no column title is sampleID\n";
    ($specimenCol >= 0) ||
	die "E parsing xlsx: no col title is specimenID\n";
    ($centerCol >= 0) ||
	die "E parsing xlsx: no column title is Center\n";
    
    foreach my $row ($rowMin+1..$rowMax) {
	my $grexome = $worksheet->get_cell($row, $grexCol)->unformatted();
	# skip "0" lines
	($grexome eq "0") && next;
	($grexome =~ /^grexome(\d\d\d\d)$/) || 
	    die "E parsing xlsx: found a grexome name that I can't parse: $grexome in row $row\n";
	my $gNum = $1;
	(defined $grex2sample[$gNum]) && 
	    die "E parsing xlsx: have 2 lines with grexomeNum $gNum\n";
	my $specimen = $worksheet->get_cell($row, $specimenCol)->value;
	my $center =  $worksheet->get_cell($row, $centerCol)->value;
	($center eq "Strasbourg") && ($specimen = "PRY-$specimen");
	$grex2sample[$gNum] = $specimen;
    }
}


#############################################
# now deal with each grexome from @grex2sample

# make sure we don't use the same infile twice (eg if globbing isn't specific enough)
# key == infile (as found by the glob), value == grexomenum if file was already used
my %infilesDone = ();

foreach my $gNum (50..$#grex2sample) {
    # silently skip grexomes that have been obsoleted (dupes)
    if (($gNum == 312) || ($gNum == 340) || ($gNum == 428) || ($gNum == 497) || ($gNum == 525) || ($gNum == 527)) {
	next;
    }
    # warn & skip if no specimen found
    (defined $grex2sample[$gNum]) || 
	((warn "W: no specimen found for grexome $gNum, skipping this grexome number\n") && next);

    my $sample = $grex2sample[$gNum];

    # precise filename patterns for each dataset:
    # strasbourg: /${sample}-R1.fastq.gz (because we prepended PRY- to $sample)
    # integragen: /Sample_${sample}_R1_fastq.gz
    # FV: /${sample}_*_R1_001.fastq.gz
    # genoscope13: /G430_CP_${sample}_[0-9]_1_*.fastq.gz
    # genoscope14 and 15: /E487_CP_${sample}_[0-9]_1_*.fastq.gz
    # novo16: /s${sample}_*_1.fq.gz
    # novo17: /P${sample}_*_1.fq.gz
    # novo19: /${sample}_*_1.fq.gz
    # Biomnis19: /${sample}_*_R1_001.fastq.gz
    # IBP_2019: /${sample}_R1_001.fastq.gz (they fucked up the names with i for 1 but I fixed it)
    # IBP_2019 second batch (01/08/2019): /${sample}_R1.fastq.gz
    # BGI_2021: /\w+_L\d_{$sample}-\d+_1.fq.gz (eg: V300096729_L4_B5EHUMazpEBAAIBAA-522_1.fq.gz)
    # these are unified into the following globs, update if needed when
    # new datasets arrive (we error out if something is fishy)
    my @files1 = glob("${inPath}/*/${sample}[-_]R1.fastq.gz ${inPath}/*/*${sample}_*R1_*.gz ${inPath}/*/*${sample}_*_1*.gz ${inPath}/*/*_${sample}-*_1.fq.gz " );
    my @files2 = glob("${inPath}/*/${sample}[-_]R2.fastq.gz ${inPath}/*/*${sample}_*R2_*.gz ${inPath}/*/*${sample}_*_2*.gz ${inPath}/*/*_${sample}-*_2.fq.gz " );

    (@files1) || 
	((warn "W: no files found for sample $sample == grexome $gNum, fix the globs, skipping this sample for now\n") && next);
    (@files1 == @files2) || 
	die "E: different numbers of files for the 2 paired-ends of sample $sample == grexome $gNum:\n@files1\n@files2\n";

    # make sure these files haven't been seen before
    foreach my $f (@files1, @files2) {
	(defined $infilesDone{$f}) &&
	    die "E: infile $f found for sample $sample == grexome $gNum but it was previously used, for grexome $infilesDone{$f}! check the logs for it! Then fix the globs\n";
	$infilesDone{$f} = $gNum;
    }

    # prepare grexome name
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";

    if ((-e "$parentDir$outDir/${grexome}_1.fq.gz") && (-e "$parentDir$outDir/${grexome}_2.fq.gz")) {
	# 09/09/2019: don't INFO when skipping, it's too much noise
	# warn "I: skipping sample $sample == grexome $grexome because targets $parentDir$outDir/${grexome}_[12].fq.gz already exist\n";
	next;
    }
    elsif ((-e "$parentDir$outDir/${grexome}_1.fq.gz") || (-e "$parentDir$outDir/${grexome}_2.fq.gz")) {
	die "E: for sample $sample == grexome $grexome , only one of the two targets $parentDir$outDir/${grexome}_[12].fq.gz already exists!\n";
    }
    # else keep going, need to make the symlinks/files



    if (@files1 == 1) {
	# prepare symlink target names
	my $f1 = $files1[0];
	my $f2 = $files2[0];
	# replace $inPath with ../$lastDir/
	($f1 =~ s~^$inPath~../$lastDir/~) || 
	    die "E: cannot replace inPath $inPath from file1 $f1 for sample $sample == grexome $gNum\n";
	($f2 =~ s~^$inPath~../$lastDir/~) || 
	    die "E: cannot replace inPath $inPath from file2 $f2 for sample $sample == grexome $gNum\n";

	my $com = "cd $parentDir$outDir/; ln -s $f1  ${grexome}_1.fq.gz ; ln -s $f2 ${grexome}_2.fq.gz" ;

	if (! $real) {
	    warn "I: dryrun, would run: $com\n";
	}
	else {
	    warn "I: single pair of files for $grexome, symlinking with: $com\n";
	    system($com);
	}
    }

    else {
	# several pairs of files, need to cat them
	my $com1 = "cat @files1 > $parentDir$outDir/${grexome}_1.fq.gz";
	my $com2 = "cat @files2 > $parentDir$outDir/${grexome}_2.fq.gz";

	if (! $real) {
	    warn "I: dryrun, would run: $com1\n";
	    warn "I: dryrun, would then run: $com2\n";
	}
	else {
	    warn "I: starting cat of $grexome with command: $com1\n";
	    system($com1);
	    warn "I: finishing cat of $grexome with command: $com2\n";
	    system($com2);
	}
    }
}

my $nbInfiles = scalar(keys(%infilesDone));
my $nbFqFiles = `cd $inPath ; /bin/ls -1 */*.gz | wc -l` ;
# some grexomes have been obsoleted because they were dupes,
# the corresponding FASTQs are still there, there are $nbObsoleteFiles,
# don't warn about them
my $nbObsoleteFiles = 12;
if ($nbInfiles + $nbObsoleteFiles == $nbFqFiles) {
    warn "\nI: examined $nbInfiles source FASTQ files and skipped $nbObsoleteFiles obsoletes in total, this is the expected number.\n";
}
else {
    warn "\nW: examined $nbInfiles source FASTQ files and skipped $nbObsoleteFiles obsoletes in total, but we actually have $nbFqFiles! why didn't we examine them all? check this!!\n";
}
