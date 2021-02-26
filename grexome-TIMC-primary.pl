#!/usr/bin/perl


# NTM
# 03/02/2021


# This is a wrapper script for the grexome-TIMC primary analysis
# pipeline, starting from "grexomized" FASTQs...
# This means that for each sample we expect a single pair of FASTQ
# files in $fastqDir, and these files must be named ${sample}_1.fq.gz
# and ${sample}_2.fq.gz .
# The "samples" are as listed in the 'sampleID' column of the provided
# $metadata file. If sampleID is '0' the row is ignored.
# Specify samples to process with --samples, otherwise every sample
# from $metadata that doesn't have a matching BAM file is processed. 
#
# Args: see $USAGE.

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Copy qw(copy move);
use File::Basename qw(basename);
use File::Temp qw(tempdir);
use File::Spec qw(splitdir);
use FindBin qw($RealBin);
use Spreadsheet::XLSX;

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


#############################################
## hard-coded paths and stuff that need to be custumized

# dir holding the hierarachy of subdirs and files containing all
# the data (FASTQs, BAMs, GVCFs). The hierarchy (specified later)
# may not need to change, but $dataDir certainly does
my $dataDir = "/data/nthierry/PierreRay/";

# number of threads / parallel jobs for fastq2bam (BWA),
# bam2gvcf* (strelka and gatk), filterBadCalls, mergeGVCFs
my $jobs = 20;

# we need 1_filterBadCalls.pl from the grexome-TIMC-Secondary repo,
# install it somewhere and set $filterBin to point to it
my $filterBin = "$RealBin/../SecondaryAnalyses/1_filterBadCalls.pl";


#############################################
## hard-coded subtrees and stuff that shouldn't need to change much

####### FASTQs
# subdir containing the "grexomized" FASTQs
my $fastqDir = "$dataDir/FASTQs_All_Grexomized/";


####### BAMs
# subdir where BAM/BAI files and associated logfiles are produced,
# this can vary depending on the run / date / server / whatever
my $bamDir = "BAMs_grexome_NTM/BAMs_NTM_Luxor/";

# subdir where all final BAMs & BAIs are symlinked
my $allBamsDir = "BAMs_All_Selected/";


####### GVCFs
# subdir where GVCF subtree is populated
my $gvcfDir = "GVCFs_grexome/";

# for each caller we produce raw GVCFs, then filter them, and finally
# merge them. Subdirs for each caller, in the order Raw - Filtered - Merged:
my @strelkaDirs = ("GVCFs_Strelka_Raw/","GVCFs_Strelka_Filtered/","GVCFs_Strelka_Filtered_Merged/");
my @gatkDirs = ("GVCFs_GATK_Raw/","GVCFs_GATK_Filtered/","GVCFs_GATK_Filtered_Merged/");
# prepend $dataDir/$gvcfDir to each
foreach my $i (0..2)  {
    $strelkaDirs[$i] = "$dataDir/$gvcfDir/".$strelkaDirs[$i];
    $gatkDirs[$i] = "$dataDir/$gvcfDir/".$gatkDirs[$i];
}


#############################################
## options / params from the command-line

# metadata file with all samples
my $metadata;

# comma-separated list of samples to process, if empty we process
# every sample from $metadata that doesn't have a BAM
my $samplesOfInterest;

# outDir must not exist, it will be created and populated
my $outDir;

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCprim_config.pm";

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Run the grexome-TIMC primary analysis pipeline, ie start from FASTQ files and:
produce BAM with fastq2bam.pl (trim, align, mark dupes, sort);
produce individual GVCFs with bam2gvcf_strelka.pl and bam2gvcf_gatk.pl;
filter low-quality variant calls with filterBadCalls.pl;
produce a merged GVCF per variant-caller with mergeGVCFs.pl.

BAMs and GVCFs are produced in a hierarchy of subdirs defined at the top of this script,
please customize them (eg \$dataDir).
Logs and copies of the metadata are produced in the provided \$outDir (which must not exist).
Each step of the pipeline is a stand-alone self-documented script, this is just a wrapper.
Every install-specific param should be in this script or in grexomeTIMCprim_config.pm.

Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--metadata string : patient metadata xlsx file, with path
--samples string : comma-separated list of sampleIDs to process, defaults to all samples in 
          metadata xlsx that don't have a BAM file
--outdir string : subdir where logs and workfiles will be created, must not pre-exist
--config string [$config] : your customized copy (with path) of the distributed *config.pm
--help : print this USAGE";

GetOptions ("metadata=s" => \$metadata,
	    "samplesOfInterest=s" => \$samplesOfInterest,
	    "outdir=s" => \$outDir,
	    "config=s" => \$config,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($metadata) || die "E $0: you must provide a metadata file\n";
(-f $metadata) || die "E $0: the supplied metadata file doesn't exist:\n$metadata\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCprim_config->import( qw(refGenome fastTmpPath) );

($outDir) || die "E $0: you must provide an outDir\n";
(-e $outDir) && 
    die "E $0: outDir $outDir already exists, remove it or choose another name.\n";


#############################################
# sanity-check all hard-coded paths (only now, so --help works)
(-d $dataDir) ||
    die "E $0: dataDir $dataDir needs to pre-exist, at least containing the FASTQs\n";

(-f $filterBin) ||
    die "E $0: cannot find filterBadCalls.pl from grexome-TIMC-Secondary, install it and set \$filterBin accordingly\n";

(-d $fastqDir) ||
    die "E $0: fastqDir $fastqDir doesn't exist\n";

(-d "$dataDir/$bamDir") || (mkdir "$dataDir/$bamDir") ||
    die "E $0: bamDir $dataDir/$bamDir doesn't exist and can't be mkdir'd\n";

(-d "$dataDir/$allBamsDir") || (mkdir "$dataDir/$allBamsDir") ||
    die "E $0: allBamsDir $dataDir/$allBamsDir doesn't exist and can't be mkdir'd\n";

(-d "$dataDir/$gvcfDir") || (mkdir "$dataDir/$gvcfDir") ||
    die "E $0: gvcfDir $dataDir/$gvcfDir doesn't exist and can't be mkdir'd\n";

foreach my $d (@strelkaDirs, @gatkDirs) {
    (-d $d) || (mkdir $d) ||
	die "E $0: a strelka or gatk GVCF subdir $d doesn't exist and can't be mkdir'd\n";
}


#########################################################
# parse patient metadata file to grab sampleIDs, limit to samples of interest

# key==existing sample of interest, value==1
my %samples = ();

{
    my $workbook = Spreadsheet::XLSX->new("$metadata");
    (defined $workbook) ||
	die "E $0: when parsing xlsx\n";
    ($workbook->worksheet_count() == 1) ||
	die "E $0: parsing xlsx: expecting a single worksheet, got ".$workbook->worksheet_count()."\n";
    my $worksheet = $workbook->worksheet(0);
    my ($colMin, $colMax) = $worksheet->col_range();
    my ($rowMin, $rowMax) = $worksheet->row_range();
    # check the column titles and grab indexes of our columns of interest
    my ($sampleCol) = (-1);
    foreach my $col ($colMin..$colMax) {
	my $cell = $worksheet->get_cell($rowMin, $col);
	# if column has no header just ignore it
	(defined $cell) || next;
	($cell->value() eq "sampleID") &&
	    ($sampleCol = $col);
    }
    ($sampleCol >= 0) ||
	die "E $0: parsing xlsx: no column title is sampleID\n";

    foreach my $row ($rowMin+1..$rowMax) {
	my $sample = $worksheet->get_cell($row, $sampleCol)->unformatted();
	# skip "0" lines
	($sample eq "0") && next;
	(defined $samples{$sample}) && 
	    die "E $0: parsing xlsx: have 2 lines with sample $sample\n";
	$samples{$sample} = 1;
    }
}

if ($samplesOfInterest) {
    # make sure every listed sample is in %samples and promote it's value to 2
    foreach my $soi (split(/,/, $samplesOfInterest)) {
	($samples{$soi}) ||
	    die "E $0: processing --samples: a specified sample $soi does not exist in the metadata file\n";
	($samples{$soi} == 1) ||
	    warn "W $0: processing --samples: sample $soi was specified twice, is that a typo?\n";
	$samples{$soi} = 2;
    }
    # now ignore all other samples
    foreach my $s (keys %samples) {
	if ($samples{$s} != 2) {
	    delete($samples{$s});
	}
    }
}

# array of sampleIDs to process, sorted (the merged GVCF will use this order)
my @samples = ();

# exclude any sample that already has a BAM, but die if called with --samples
foreach my $s (sort keys(%samples)) {
    my $bam = "$dataDir/$allBamsDir/$s.bam";
    if (! -e $bam) {
	push(@samples, $s);
    }
    elsif ($samples{$s} == 2) {
	die "E $0: sample $s specified with --samples already has a BAM: $bam\n";
    }
}

# comma-separated string of samples to process
my $samples = join(',',@samples);
    

#############################################

my $now = strftime("%F %T", localtime);
warn "I $0: $now - starting to run\n\n";


# prep is AOK, we can mkdir outDir now
mkdir($outDir) || die "E $0: cannot mkdir outDir $outDir\n";

# copy the provided metadata file into $outDir
copy($metadata, $outDir) ||
    die "E $0: cannot copy metadata to outDir: $!\n";
# use the copied versions in scripts (eg if original gets edited while analysis is running...)
$metadata = "$outDir/".basename($metadata);


# randomly-named subdir of &fastTmpPath() (to avoid clashes),
# $tmpDir is removed afterwards
my $tmpDir = tempdir(DIR => &fastTmpPath());

################################
# BAMS

# make BAMs
my $com = "perl $RealBin/2_Fastq2Bam/fastq2bam.pl --indir $fastqDir --samples $samples --outdir $dataDir/$bamDir ";
$com .= "--genome ".&refGenome()." --threads $jobs --real ";
system($com) && die "E $0: fastq2bam FAILED: $!";
$now = strftime("%F %T", localtime);
warn "I $0: $now - fastq2bam DONE, stepwise logfiles are available as $dataDir/$bamDir/*log\n\n";

# symlink just the BAMs/BAIs in $allBamsDir with relative symlinks (so rsyncing the
# whole tree elsewhere still works)
# building the relative path corrctly is a bit tricky
{
    my @bDirs = splitdir($bamDir);
    my @abDirs = splitdir($allBamsDir);
    # remove common leading dirs
    while($bDirs[0] eq $abDirs[0]) {
	shift(@bDirs);
	shift(@abDirs);
    }
    # build relative path from allBamsDir to bamDir
    my $relPath = '../' x scalar(@abDirs);
    $relPath .= join('/',@bDirs);
    foreach my $sample (@samples) {
	foreach my $file ("$sample.bam", "$sample.bam.bai") {
	    (! -f "$dataDir/$bamDir/$file") &&
		(warn "W $0: I want to symlink $dataDir/$bamDir/$file but it doesn't exist, skipping\n") && next;
	    (-e "$dataDir/$allBamsDir/$file") &&
		(warn "W $0: I want to symlink $dataDir/$bamDir/$file but symlink already exists, skipping\n") && next;
	    symlink("$relPath/$file", "$dataDir/$allBamsDir/$file") ||
		die "E $0: cannot symlink $relPath/$file : $!";
	}
    }
}
$now = strftime("%F %T", localtime);
warn "I $0: $now - symlinking BAMs/BAIs in $allBamsDir DONE\n\n";


################################
# STRELKA

# make INDIVIDUAL STRELKA GVCFs
my $strelkaWorkdir =  "$outDir/StrelkaResults/";
$com = "perl $RealBin/3_Bam2Gvcf_Strelka/bam2gvcf_strelka.pl --indir $dataDir/$allBamsDir --samples $samples --outdir $strelkaWorkdir --jobs $jobs --config $config --real";
system($com) && die "E $0: bam2gvcf_strelka FAILED: $?";
# check logs:
open(LOGS, "cat $strelkaWorkdir/*/workflow.exitcode.txt |") ||
    die "E $0: cannot open strelka logs: $!\n";
while (my $line = <LOGS>) {
    chomp($line);
    ($line eq "0") ||
	die "E $0: non-zero exit code from a strelka run, check $strelkaWorkdir/*/workflow.exitcode.txt\n";
}
close(LOGS);

# move STRELKA GVCFs and TBIs into $gvcfDir subtree
$com = "perl $RealBin/3_Bam2Gvcf_Strelka/moveGvcfs.pl $strelkaWorkdir $strelkaDirs[0]";
system($com) && die "E $0: strelka moveGvcfs FAILED: $?";

$now = strftime("%F %T", localtime);
warn "I $0: $now - variant-calling with STRELKA DONE\n\n";


# filter INDIVIDUAL STRELKA GVCFs
foreach my $sample (@samples) {
    warn "I $0: starting filter of strelka $sample\n";
    $com = "bgzip -cd -\@6 $strelkaDirs[0]/$sample.g.vcf.gz | ";
    $com .= "perl $filterBin --metadata=$metadata --tmpdir=$tmpDir/Filter --keepHR --jobs $jobs | ";
    $com .= "bgzip -c -\@2 > $strelkaDirs[1]/$sample.filtered.g.vcf.gz";
    system($com) && die "E $0: filterGVCFs for strelka $sample FAILED: $?";
}
$now = strftime("%F %T", localtime);
warn "I $0: $now - filtering strelka GVCFs DONE\n\n";


################################
# GATK

# make INDIVIDUAL GATK GVCFs
my $gatkWorkdir = "$outDir/GatkResults/";
$com = "perl $RealBin/3_Bam2Gvcf_GATK/bam2gvcf_gatk.pl --indir $dataDir/$allBamsDir --samples $samples --outdir $gatkWorkdir --jobs $jobs --config $config --real";
system($com) && die "E $0: bam2gvcf_gatk FAILED: $?";

# move GATK GVCFs + TBIs + logs into $gvcfDir subtree and remove now-empty workdir:
foreach my $sample (@samples) {
    foreach my $file ("$sample.g.vcf.gz", "$sample.g.vcf.gz.tbi", "$sample.log") {
	move("$gatkWorkdir/$file", "$gatkDirs[0]") ||
	    die "E $0: cannot move $gatkWorkdir/$file to $gatkDirs[0] : $!";
    }
}
rmdir($gatkWorkdir) || die "E $0: cannot rmdir gatkWorkdir $gatkWorkdir: $!";

$now = strftime("%F %T", localtime);
warn "I $0: $now - variant-calling with GATK DONE\n\n";

# filter INDIVIDUAL GATK GVCFs
foreach my $sample (@samples) {
    warn "I $0: starting filter of gatk $sample\n";
    $com = "bgzip -cd -\@6 $gatkDirs[0]/$sample.g.vcf.gz | ";
    $com .= "perl $filterBin --metadata=$metadata --tmpdir=$tmpDir/Filter --keepHR --jobs $jobs | ";
    $com .= "bgzip -c -\@2 > $gatkDirs[1]/$sample.filtered.g.vcf.gz";
    system($com) && die "E $0: filterGVCFs for gatk $sample FAILED: $?";
}
$now = strftime("%F %T", localtime);
warn "I $0: $now - filtering gatk GVCFs DONE\n\n";


################################
# merge new GVCFs with previous merged

# YYMMDD for creating timestamped new merged
my $date = strftime("%y%m%d", localtime);


# make batchfile with list of GVCFs to merge
my $batchFile = "$outDir/batchFile_strelka.txt";
open(BATCH, ">$batchFile") ||
    die "E $0: cannot create strelka batchFile $batchFile: $!\n";

# we want to merge the new GVCFs with the most recent previous merged,
# a bit ugly but functional
my $prevMerged = `ls -rt1 "$strelkaDirs[2]/*.g.vcf.gz" | tail -n 1`;
print BATCH "$prevMerged\n";
foreach my $sample (@samples) {
    print BATCH "$strelkaDirs[1]/$sample.filtered.g.vcf.gz\n";
}
close(BATCH);

my $newMerged = "$strelkaDirs[2]/grexomes_strelka_merged_$date.g.vcf.gz";
(-e $newMerged) &&
    die "E $0: want to merge GVCFs but newMerged already exists: $newMerged\n";

# -> merge:
$com = "perl $RealBin/4_MergeGVCFs/mergeGVCFs.pl --filelist $batchFile --config $config --cleanheaders --jobs $jobs ";
$com .= "2>  $outDir/merge_strelka.log ";
$com .= "| bgzip -c -\@12 > $newMerged";
warn "I $0: starting to merge strelka GVCFs\n";
system($com) && die "E $0: mergeGvcfs for strelka FAILED: $?";
$now = strftime("%F %T", localtime);
warn "I $0: $now - merging strelka GVCFs DONE\n\n";

# index
$com = "tabix $newMerged";
system($com) && die "E $0: tabix for merged strelka FAILED: $?";
$now = strftime("%F %T", localtime);
warn "I $0: $now - indexing merged strelka GVCF DONE\n\n";


# same for GATK: merge new GVCFs with most recent previous merged
$batchFile = "$outDir/batchFile_gatk.txt";
open(BATCH, ">$batchFile") ||
    die "E $0: cannot create gatk batchFile $batchFile: $!\n";
$prevMerged = `ls -rt1 "$gatkDirs[2]/*.g.vcf.gz" | tail -n 1`;
print BATCH "$prevMerged\n";
foreach my $sample (@samples) {
    print BATCH "$gatkDirs[1]/$sample.filtered.g.vcf.gz\n";
}
close(BATCH);

$newMerged = "$gatkDirs[2]/grexomes_gatk_merged_$date.g.vcf.gz";
(-e $newMerged) &&
    die "E $0: want to merge GVCFs but newMerged already exists: $newMerged\n";

# -> merge:
$com = "perl $RealBin/4_MergeGVCFs/mergeGVCFs.pl --filelist $batchFile --config $config --cleanheaders --jobs $jobs ";
$com .= "2>  $outDir/merge_gatk.log ";
$com .= "| bgzip -c -\@12 > $newMerged";
warn "I $0: starting to merge gatk GVCFs\n";
system($com) && die "E $0: mergeGvcfs for gatk FAILED: $?";
$now = strftime("%F %T", localtime);
warn "I $0: $now - merging gatk GVCFs DONE\n\n";

# index
$com = "tabix $newMerged";
system($com) && die "E $0: tabix for merged gatk FAILED: $?";
$now = strftime("%F %T", localtime);
warn "I $0: $now - indexing merged gatk GVCF DONE\n\n";


warn "I $0: ALL DONE, please examine the logs and if AOK you can remove\n";
warn "I $0: obsolete merged GVCFs and sync all results to cargo:bettik\n";
warn "I $0: with the following commands:\n";

my $oldestMergedStrelka = `ls -rt1 "$strelkaDirs[2]/*.g.vcf.gz" | head -n 1`;
my $oldestMergedGatk = `ls -rt1 "$gatkDirs[2]/*.g.vcf.gz" | head -n 1`;

$com = "cd $dataDir\n";
$com .= "rm -i $oldestMergedStrelka $oldestMergedStrelka.tbi\n";
$com .= "rm -i $oldestMergedGatk $oldestMergedGatk.tbi\n";
$com .= "rsync -rtvn --delete $bamDir/ cargo:/bettik/thierryn/$bamDir/\n";
$com .= "rsync -rtvln --delete $allBamsDir/ cargo:/bettik/thierryn/$allBamsDir/\n";
$com .= "rsync -rtvn --delete $gvcfDir/ cargo:/bettik/thierryn/$gvcfDir/\n";
$com .= "## redo without -n if AOK:\n";
$com .= "rsync -rtv --delete $bamDir/ cargo:/bettik/thierryn/$bamDir/\n";
$com .= "rsync -rtvl --delete $allBamsDir/ cargo:/bettik/thierryn/$allBamsDir/\n";
$com .= "rsync -rtv --delete $gvcfDir/ cargo:/bettik/thierryn/$gvcfDir/\n";

warn "$com\n";
