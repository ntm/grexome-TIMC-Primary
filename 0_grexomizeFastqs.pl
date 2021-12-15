#!/usr/bin/perl

# 13/06/2019
# NTM

# This is a house-keeping script used to organize our FASTQ files before
# running the grexome-TIMC-primary pipeline. It's probably not re-usable as-is.
#
# Sometimes we get several pairs of FASTQ files for a single sample,
# eg when the sample was sequenced on several lanes.
# Also, the naming conventions are always different.
#
# This script creates a single pair of "files" (symlinks when possible)
# in $inPath/../$outDir/ for each sample, with uniform naming scheme
# ${sample}_[12].fq.gz .
#
# Basically we want to symlink when there's a single pair of files,
# and otherwise cat @infiles1 > $outDir/${sample}_1.fq.gz
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
use File::Basename qw(basename);
use FindBin qw($RealBin);
use Getopt::Long;

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);

# make glob case-insensitive (for eg */P17if085*)
use File::Glob qw(:globally :nocase);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# filename with path to samples metadata xlsx file
my $samplesFile;

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

my $USAGE = "
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile string : samples metadata xlsx file, with path.
--inpath string : path to a directory containing subdirs, each containing the
    actual original fastq.gz files (eg FASTQs_grexome/).
--outdir string [default: $outDir] : name of brother-dir of \$inPath where \$sampleID_*.fq.gz
    files are symlinked/produced. NOTE: outDir cannot contain a slash, it will be created
    if needed and will be a BROTHER of the final subdir of \$inPath ie \$inPath/../\$outDir.
--real : actually do the work, default is a dry run ie only print info on what would be done.
--help : print this USAGE";


GetOptions ("samplesFile=s" => \$samplesFile,
	    "inpath=s" => \$inPath,
	    "outdir=s" => \$outDir,
	    "real" => \$real,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n\n$USAGE\n");


# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samplesFile file\n";
(-f $samplesFile) || die "E $0: the supplied samplesFile file doesn't exist\n";

($inPath) || 
    die "E: you MUST specify the path holding subdirs that contain the original FASTQs, with --inpath\n";
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
## parse $samplesFile

# $sample2specimenR : hashref, key is sampleID, value is specimenID
my $sample2specimenR;

{
    my @parsed = &parseSamples($samplesFile);
    $sample2specimenR = $parsed[1];
}


#############################################
# now deal with each sample

# make sure we don't use the same infile twice (eg if globbing isn't specific enough)
# key == infile (as found by the glob), value == $sample if file was already used
my %infilesDone = ();

foreach my $sample (sort keys(%$sample2specimenR)) {
    my $specimen = $sample2specimenR->{$sample};

    # precise filename patterns for each dataset:
    # strasbourg: /${specimen}-R1.fastq.gz
    # integragen: /Sample_${specimen}_R1_fastq.gz
    # FV: /${specimen}_*_R1_001.fastq.gz
    # genoscope13: /G430_CP_${specimen}_[0-9]_1_*.fastq.gz
    # genoscope14 and 15: /E487_CP_${specimen}_[0-9]_1_*.fastq.gz
    # novo16: /s${specimen}_*_1.fq.gz
    # novo17: /P${specimen}_*_1.fq.gz
    # novo19: /${specimen}_*_1.fq.gz
    # Biomnis19: /${specimen}_*_R1_001.fastq.gz
    # IBP_2019: /${specimen}_R1_001.fastq.gz (they fucked up the names with i for 1 but I fixed it)
    # IBP_2019 second batch (01/08/2019): /${specimen}_R1.fastq.gz
    # BGI_2021: /\w+_L\d_{$specimen}-\d+_1.fq.gz (eg: V300096729_L4_B5EHUMazpEBAAIBAA-522_1.fq.gz)
    # Macrogen_2021: /${specimen}_1.fastq.gz
    # these are unified into the following globs, update if needed when
    # new datasets arrive (we error out if something is fishy)
    my @files1 = glob("${inPath}/*/${specimen}[-_]R1.fastq.gz ${inPath}/*/*${specimen}_*R1_*.gz ${inPath}/*/*${specimen}_*_1*.gz ${inPath}/*/*_${specimen}-*_1.fq.gz ${inPath}/*/${specimen}_1.fastq.gz " );
    my @files2 = glob("${inPath}/*/${specimen}[-_]R2.fastq.gz ${inPath}/*/*${specimen}_*R2_*.gz ${inPath}/*/*${specimen}_*_2*.gz ${inPath}/*/*_${specimen}-*_2.fq.gz ${inPath}/*/${specimen}_2.fastq.gz " );

    (@files1) || 
	((warn "W: no files found for $sample == specimen $specimen, fix the globs, skipping this sample for now\n") && next);
    (@files1 == @files2) || 
	die "E: different numbers of files for the 2 paired-ends of $sample == specimen $specimen:\n@files1\n@files2\n";

    # make sure these files haven't been seen before
    foreach my $f (@files1, @files2) {
	(defined $infilesDone{$f}) &&
	    die "E: infile $f found for $sample == specimen $specimen, but it was previously used for $infilesDone{$f}! check the logs for it! Then fix the globs\n";
	$infilesDone{$f} = $sample;
    }

    if ((-e "$parentDir$outDir/${sample}_1.fq.gz") && (-e "$parentDir$outDir/${sample}_2.fq.gz")) {
	# 09/09/2019: don't INFO when skipping, it's too much noise
	# warn "I: skipping $sample == specimen $specimen because grexomized fastqs already exist\n";
	next;
    }
    elsif ((-e "$parentDir$outDir/${sample}_1.fq.gz") || (-e "$parentDir$outDir/${sample}_2.fq.gz")) {
	die "E: $sample already has one grexomized FASTQ but not the other! Something is wrong, investigate!\n";
    }
    # else keep going, need to make the symlinks/files

    if (@files1 == 1) {
	# prepare symlink target names
	my $f1 = $files1[0];
	my $f2 = $files2[0];
	# replace $inPath with ../$lastDir/
	($f1 =~ s~^$inPath~../$lastDir/~) || 
	    die "E: cannot replace inPath $inPath from file1 $f1 for specimen $specimen == $sample\n";
	($f2 =~ s~^$inPath~../$lastDir/~) || 
	    die "E: cannot replace inPath $inPath from file2 $f2 for specimen $specimen == $sample\n";

	my $com = "cd $parentDir$outDir/; ln -s $f1 ${sample}_1.fq.gz ; ln -s $f2 ${sample}_2.fq.gz" ;

	if (! $real) {
	    warn "I: dryrun, would run: $com\n";
	}
	else {
	    warn "I: single pair of files for $sample, symlinking with: $com\n";
	    system($com);
	}
    }

    else {
	# several pairs of files, need to cat them
	my $com1 = "cat @files1 > $parentDir$outDir/${sample}_1.fq.gz";
	my $com2 = "cat @files2 > $parentDir$outDir/${sample}_2.fq.gz";

	if (! $real) {
	    warn "I: dryrun, would run: $com1\n";
	    warn "I: dryrun, would then run: $com2\n";
	}
	else {
	    warn "I: starting cat of $sample with command: $com1\n";
	    system($com1);
	    warn "I: finishing cat of $sample with command: $com2\n";
	    system($com2);
	}
    }
}

my $nbInfiles = scalar(keys(%infilesDone));
my $nbFqFiles = `cd $inPath ; /bin/ls -1 */*.gz | wc -l` ;
chomp($nbFqFiles);
# some grexomes have been obsoleted because they were dupes,
# the corresponding FASTQs are still there, don't warn about them.
# number of obsolete FASTQ files: hard-coded here, 
my $nbObsoleteFiles = 12;
if ($nbInfiles + $nbObsoleteFiles == $nbFqFiles) {
    warn "\nI: nb of examined ($nbInfiles) + skipped obsolete ($nbObsoleteFiles) FASTQ files == nbFqFiles found with ls|wc, good!\n";
}
else {
    warn "\nW: examined $nbInfiles FASTQs and skipped $nbObsoleteFiles obsoletes, but we actually have $nbFqFiles! why didn't we examine them all? check this!!\n";
}
