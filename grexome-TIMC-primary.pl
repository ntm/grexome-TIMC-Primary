#!/usr/bin/perl


# NTM
# 03/02/2021


# This is a wrapper script for the grexome-TIMC primary analysis
# pipeline, starting from "grexomized" FASTQs...
# This means that for each sample we expect a single pair of FASTQ
# files in $fastqDir, and these files must be named ${sample}_1.fq.gz
# and ${sample}_2.fq.gz .
# The "samples" are as listed in the 'sampleID' column of the provided
# $samplesFile. If sampleID is '0' the row is ignored.
# You can specify samples to process with --samplesOfInterest, otherwise 
# every sample from $samplesFile is processed. 
#
# Args: see $USAGE.

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Copy qw(copy move);
use File::Path qw(remove_tree);
use File::Basename qw(basename);
use File::Temp qw(tempdir);
use File::Spec;
use FindBin qw($RealBin);

use lib "$RealBin";
use grexome_metaParse qw(parseSamples);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## hard-coded subtrees and stuff that shouldn't need to change much

####### BAMs
# $bamDir: subdir of $dataDir where BAM/BAI files and associated logfiles
# are produced, this can vary depending on the run / date / server / whatever...
my $bamDir = "BAMs_grexome/";
# ..., but in any case we always symlink all final BAMs & BAIs in $dataDir/$allBamsDir
my $allBamsDir = "BAMs_All_Selected/";


####### GVCFs
# subdir of $dataDir where GVCF subtree is populated
my $gvcfDir = "GVCFs_grexome/";

# for each caller we produce raw GVCFs, then filter them, and finally
# merge them.
# $callerDirs{$caller} is a ref to an array of 3 dirs for $caller, in the
# order Raw - Filtered - Merged ($dataDir/$gvcfDir will be prepended soon):
my %callerDirs = (
    "strelka" => ["GVCFs_Strelka_Raw/","GVCFs_Strelka_Filtered/","GVCFs_Strelka_Filtered_Merged/"],
    "gatk" => ["GVCFs_GATK_Raw/","GVCFs_GATK_Filtered/","GVCFs_GATK_Filtered_Merged/"]);


#############################################
## programs that we use
my $bgzip = "bgzip";
my $tabix = "tabix";
my $zgrep = "zgrep";
(`which $bgzip` =~ /$bgzip$/) ||
    die "E $0: the bgzip executable $bgzip can't be found\n";
(`which $tabix` =~ /$tabix$/) ||
    die "E $0: the tabix executable $tabix can't be found\n";
(`which $zgrep` =~ /$zgrep$/) ||
    die "E $0: the zgrep executable $zgrep can't be found\n";


#############################################
## options / params from the command-line

# samples metada file
my $samplesFile;

# comma-separated list of samples of interest, if empty we process
# every sample from $samplesFile (skipping any step where the resulting
# outfile already exists)
my $SOIs;

# workDir must not exist, it will be created and populated
my $workDir;

# indication of the number of threads / parallel jobs that we can run
my $jobs = 20;

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCprim_config.pm";

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Run the grexome-TIMC primary analysis pipeline, ie start from FASTQ files and:
- produce BAMs with fastq2bam.pl (trim, align, mark dupes, sort);
- for each specified variant-caller:
     produce individual GVCFs with bam2gvcf_\$caller.pl;
     filter low-quality variant calls with filterBadCalls.pl;
     produce a merged GVCF per variant-caller with mergeGVCFs.pl.

BAMs and GVCFs are produced in a hierarchy of subdirs defined at the top of this script,
all within \$dataDir which you must customize. The hierarchy can also be modified if desired.
Logs and copies of the metadata are produced in the provided \$workDir (which must not exist).
Each step of the pipeline is a stand-alone self-documented script, this is just a wrapper.
For each sample, any step where the result file already exists is skipped.
Every install-specific param should be in this script or in grexomeTIMCprim_config.pm.

Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile string : samples metadata xlsx file, with path
--samplesOfInterest string : comma-separated list of sampleIDs to process, default = all samples in samplesFile
--workdir string : subdir where logs and workfiles will be created, must not pre-exist
--jobs int [$jobs] : approximate number of threads/jobs that we should run
--config string [$config] : your customized copy (with path) of the distributed *config.pm
--help : print this USAGE";

GetOptions ("samplesFile=s" => \$samplesFile,
	    "samplesOfInterest=s" => \$SOIs,
	    "workdir=s" => \$workDir,
	    "jobs=i" => \$jobs,
	    "config=s" => \$config,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samples metadata file. Try $0 --help\n";
(-f $samplesFile) || die "E $0: the supplied samples metadata file doesn't exist: $samplesFile\n";

# immediately import $config, so we die if file is broken
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCprim_config->import( qw(dataDir fastqDir mirror refGenome fastTmpPath) );

($workDir) || die "E $0: you must provide a workDir. Try $0 --help\n";
(-e $workDir) && 
    die "E $0: workDir $workDir already exists, remove it or choose another name.\n";


#############################################
# sanity-check all hard-coded paths (only now, so --help works)

# dir holding the hierarachy of subdirs and files containing all
# the data (FASTQs, BAMs, GVCFs). The hierarchy (specified later)
# may not need to change, but $dataDir certainly does
my $dataDir = &dataDir();

# dir containing the "grexomized" FASTQs
my $fastqDir = &fastqDir();

(-d $dataDir) || (mkdir "$dataDir") ||
    die "E $0: dataDir $dataDir doesn't exist and can't be mkdir'd\n";

(-d $fastqDir) ||
    die "E $0: fastqDir $fastqDir doesn't exist\n";

(-d "$dataDir/$bamDir") || (mkdir "$dataDir/$bamDir") ||
    die "E $0: bamDir $dataDir/$bamDir doesn't exist and can't be mkdir'd\n";

(-d "$dataDir/$allBamsDir") || (mkdir "$dataDir/$allBamsDir") ||
    die "E $0: allBamsDir $dataDir/$allBamsDir doesn't exist and can't be mkdir'd\n";

(-d "$dataDir/$gvcfDir") || (mkdir "$dataDir/$gvcfDir") ||
    die "E $0: gvcfDir $dataDir/$gvcfDir doesn't exist and can't be mkdir'd\n";

# now that $dataDir exists, prepend $dataDir/$gvcfDir to %callerDirs and sanity-check
foreach my $caller (keys %callerDirs) {
    foreach my $i (0..2) {
	my $d = "$dataDir/$gvcfDir/".$callerDirs{$caller}->[$i];
	(-d $d) || (mkdir $d) ||
	    die "E $0: GVCF subdir $d for caller $caller doesn't exist and can't be mkdir'd\n";
	$callerDirs{$caller}->[$i] = $d;
    }
}

#########################################################
# parse samples metadata file to grab sampleIDs, limit to --samplesOfInterest if specified

# key==existing sample to process
my %samples = ();

# just use the first hashref from parseSamples, ignoring pathologyIDs
my ($s2pathoR) = &parseSamples($samplesFile);
foreach my $s (keys %$s2pathoR) {
    $samples{$s} = 1;
}

if ($SOIs) {
    # make sure every listed sample is in %samples and promote it's value to 2
    foreach my $soi (split(/,/, $SOIs)) {
	($samples{$soi}) ||
	    die "E $0: processing --samplesOfInterest: a specified sample $soi does not exist in the samples metadata file\n";
	($samples{$soi} == 1) ||
	    warn "W $0: processing --samplesOfInterest: sample $soi was specified twice, is that a typo? Skipping the dupe\n";
	$samples{$soi} = 2;
    }
    # now ignore all other samples
    foreach my $s (keys %samples) {
	if ($samples{$s} != 2) {
	    delete($samples{$s});
	}
    }
}

# exclude any sample that doesn't have FASTQs, but die if called with --samplesOfInterest
foreach my $s (sort(keys %samples)) {
    my $f1 = "$fastqDir/${s}_1.fq.gz";
    my $f2 = "$fastqDir/${s}_2.fq.gz";
    if ((! -f $f1) || (! -f $f2)) {
	if ($samples{$s} == 2) {
	    die "E $0: sample $s from --samplesOfInterest doesn't have FASTQs (looking for $f1 and $f2)\n";
	}
	else {
	    warn "W $0: sample $s from samplesFile doesn't have FASTQs, skipping it\n";
	    delete($samples{$s});
	}
    }
}

#############################################

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - starting to run\n";


# prep is AOK, we can mkdir workDir now
mkdir($workDir) || die "E $0: cannot mkdir workDir $workDir\n";

# copy the provided samples metadata file into $workDir
copy($samplesFile, $workDir) ||
    die "E $0: cannot copy metadata to workDir: $!\n";
# use the copied versions in scripts (eg if original gets edited while analysis is running...)
$samplesFile = "$workDir/".basename($samplesFile);


# randomly-named subdir of &fastTmpPath() (to avoid clashes),
# $tmpDir is removed afterwards
my $tmpDir = tempdir(DIR => &fastTmpPath(), CLEANUP => 1);

################################
# MAKE BAMS

# samples to process: those without BAMs in $dataDir/$bamDir
my $samples = "";
foreach my $s (sort(keys %samples)) {
    my $bam = "$dataDir/$bamDir/$s.bam";
    (-e $bam) || ($samples .= "$s,");
}
if ($samples) {
    # remove trailing ','
    (chop($samples) eq ',') ||
	die "E $0 chopped samples isn't ',' impossible\n";
    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - fastq2bam will process sample(s) $samples\n";
    # make BAMs
    my $com = "perl $RealBin/1_fastq2bam.pl --indir $fastqDir --samples $samples --outdir $dataDir/$bamDir ";
    $com .= "--genome ".&refGenome()." --threads $jobs --real ";
    system($com) && die "E $0: fastq2bam FAILED: $!";
    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - fastq2bam DONE, logfiles are available as $dataDir/$bamDir/*log\n";
}
else {
    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - fastq2bam DONE, all samples already had BAMs\n";
}

################################
# SYMLINK BAMS

# symlink the BAMs/BAIs in $allBamsDir with relative symlinks (so rsyncing the
# whole tree elsewhere still works)
# samples to process: those without BAMs in $dataDir/$allBamsDir
$samples = "";
foreach my $s (sort(keys %samples)) {
    my $bam = "$dataDir/$allBamsDir/$s.bam";
    (-e $bam) || ($samples .= "$s,");
}
if ($samples) {
    # remove trailing ','
    (chop($samples) eq ',') ||
	die "E $0 chopped samples isn't ',' impossible\n";
   
    # building the relative path correctly is a bit tricky
    {
	my @bDirs = File::Spec->splitdir($bamDir);
	my @abDirs = File::Spec->splitdir($allBamsDir);
	# remove common leading dirs
	while($bDirs[0] eq $abDirs[0]) {
	    shift(@bDirs);
	    shift(@abDirs);
	}
	# remove last dir if empty  (happens if eg $bamDir was slash-terminated)
	($bDirs[$#bDirs]) || pop(@bDirs);
	($abDirs[$#abDirs]) || pop(@abDirs);
	# build relative path from allBamsDir to bamDir
	my $relPath = '../' x scalar(@abDirs);
	$relPath .= join('/',@bDirs);
	foreach my $s (split(/,/,$samples)) {
	    foreach my $file ("$s.bam", "$s.bam.bai") {
		(-e "$dataDir/$bamDir/$file") ||
		    die "E $0: BAM/BAI doesn't exist but should have been made: $dataDir/$bamDir/$file\n";
		symlink("$relPath/$file", "$dataDir/$allBamsDir/$file") ||
		    die "E $0: cannot symlink $relPath/$file : $!";
	    }
	}
    }
}
$now = strftime("%F %T", localtime);
warn "I $now: $0 - symlinking any new BAMs/BAIs in $allBamsDir DONE\n";


################################
# GVCFs

# YYMMDD for creating timestamped new merged GVCFs
# set now so all callers use the same timestamp even if gatk takes a lot longer
my $date = strftime("%y%m%d", localtime);

# PIDs from tabix-index subprocesses
my @childrenPids;

# mostly same code for all callers, except house-keeping
foreach my $caller (sort(keys %callerDirs)) {
    ################################
    # make INDIVIDUAL GVCFs
    # samples to process: those without a raw GVCF
    $samples = "";
    foreach my $s (sort(keys %samples)) {
	my $gvcf = $callerDirs{$caller}->[0]."/$s.g.vcf.gz";
	(-e $gvcf) || ($samples .= "$s,");
    }
    if ($samples) {
	# remove trailing ','
	(chop($samples) eq ',') ||
	    die "E $0 chopped samples isn't ',' impossible\n";

	my $callerWorkDir =  "$workDir/Results_$caller/";

	# each caller-specific wrapper script MUST BE named precisely like this:
	my $b2gBin = "$RealBin/2_bam2gvcf_$caller.pl";
	(-e $b2gBin) ||
	    die "E $0: trying to bam2gvcf for $caller, but  b2gBin $b2gBin doesn't exist\n";
	my $com = "perl $b2gBin --indir $dataDir/$allBamsDir --samples $samples --outdir $callerWorkDir --jobs $jobs --config $config --real";
	system($com) && die "E $0: bam2gvcf_$caller FAILED: $?";

	##################
	# caller-specific code: log-checking and house-keeping
	if ($caller eq "strelka") {
	    # check that logs are empty
	    foreach my $s (split(/,/,$samples)) {
		if (! -z "$callerWorkDir/$s/workflow.error.log.txt") {
		    die "E $0: non-empty strelka error log for $s, go look in $callerWorkDir/$s/\n";
		}
		elsif (! -z "$callerWorkDir/$s/workflow.warning.log.txt") {
		    warn"W $0: non-empty strelka warning log for $s, go look in $callerWorkDir/$s/\n";
		}
		else {
		    # remove useless strelka left-overs
		    remove_tree("$callerWorkDir/$s/workspace/pyflow.data/logs/tmp/");
		}
	    }
	    # move STRELKA GVCFs and TBIs into $gvcfDir subtree
	    $com = "perl $RealBin/2_bam2gvcf_strelka_moveGvcfs.pl $callerWorkDir ".$callerDirs{"strelka"}->[0];
	    system($com) && die "E $0: strelka moveGvcfs FAILED: $?";
	}
	elsif ($caller eq "gatk") {
	    # GATK logs are a mess: they seem to adopt a format but then don't respect it,
	    # not even trying to parse it
	    # move GATK GVCFs + TBIs + logs into $gvcfDir subtree and remove now-empty callerWorkDir:
	    foreach my $s (split(/,/,$samples)) {
		foreach my $file ("$s.g.vcf.gz", "$s.g.vcf.gz.tbi", "$s.log") {
		    move("$callerWorkDir/$file", $callerDirs{"gatk"}->[0]) ||
			die "E $0: cannot move $callerWorkDir/$file to ".$callerDirs{"gatk"}->[0]." : $!";
		}
	    }
	    rmdir($callerWorkDir) || die "E $0: cannot rmdir gatk callerWorkDir $callerWorkDir: $!";
	}
	else {
	    die "E $0: new caller $caller, need to implement log-checking and house-keeping after bam2gvcf\n";
	}
	# end of caller-specific code
	##################
	
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - variant-calling with $caller DONE\n";
    }
    else {
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - variant-calling with $caller DONE, all samples already had raw $caller GVCFs\n";
    }

    ################################
    # filter INDIVIDUAL GVCFs
    foreach my $s (sort(keys %samples)) {
	# samples to process: those without a filtered GVCF
	my $gvcf = $callerDirs{$caller}->[1]."/$s.filtered.g.vcf.gz";
	(-e $gvcf) && next;
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - starting filter of $caller $s\n";
	my $com = "$bgzip -cd -\@6 ".$callerDirs{$caller}->[0]."/$s.g.vcf.gz | ";
	# giving 6+2 threads to bgzip, so reduce $jobs for filterBin: max(j/2, j-4)
	my ($jf1,$jf2) = (int(($jobs+1)/2), $jobs-4);
	my $jobsFilter = $jf1;
	($jf2 > $jf1) && ($jobsFilter = $jf2);
	$com .= "perl $RealBin/3_filterBadCalls.pl --samplesFile=$samplesFile --tmpdir=$tmpDir/Filter --keepHR --jobs $jobsFilter | ";
	$com .= "$bgzip -c -\@2 > $gvcf";
	system($com) && die "E $0: filterGVCFs for $caller $s FAILED: $?";
    }
    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - filtering $caller GVCFs DONE\n";

    ################################
    # merge new GVCFs with the most recent previous merged if one existed

    # samples in prevMerged, to avoid dupes
    my %samplesPrev;

    # make batchfile with list of GVCFs to merge
    my $batchFile = "$workDir/batchFile_$caller.txt";
    open(BATCH, ">$batchFile") ||
	die "E $0: cannot create $caller batchFile $batchFile: $!\n";

    # we want to merge the new GVCFs with the most recent previous merged,
    # if there was one. code is a bit ugly but functional
    my $prevMerged = `ls -rt1 $callerDirs{$caller}->[2]/*.g.vcf.gz | tail -n 1`;
    chomp($prevMerged);
    if ($prevMerged) {
	open(CHR, "$zgrep -m 1 ^#CHROM $prevMerged |") ||
	    die "E $0: cannot zgrep #CHROM line in prevMerged $prevMerged\n";
	my $header = <CHR>;
	chomp($header);
	my @header = split(/\t/,$header);
	# first 9 fields are regular VCF headers CHROM->FORMAT
	splice(@header, 0, 9);
	foreach my $s (@header) {
	    $samplesPrev{$s} &&
		die "E $0: most recent merged $caller GVCF has dupe sample $s! Investigate $prevMerged\n";
	    $samplesPrev{$s} = 1;
	}
	close(CHR);
	print BATCH "$prevMerged\n";
    }

    # only merge if there's at least one new sample
    my $doMerge = 0;
    foreach my $s (sort(keys %samples)) {
	# only merge $s if it's not already in prevMerged
	if (! $samplesPrev{$s}) {
	    print BATCH $callerDirs{$caller}->[1]."/$s.filtered.g.vcf.gz\n";
	    $doMerge = 1;
	}
    }
    close(BATCH);

    if ($doMerge) {
	my $newMerged = $callerDirs{$caller}->[2]."/grexomes_${caller}_merged_$date.g.vcf.gz";
	(-e $newMerged) &&
	    die "E $0: want to merge GVCFs but newMerged already exists: $newMerged\n";

	# -> merge:
	# NOTE we are piping to bgzip -@12 , so we exceed $jobs quite a bit.
	# if this turns into a problem we can tune it down (eg $jobsFilter)
	my $com = "perl $RealBin/4_mergeGVCFs.pl --filelist $batchFile --config $config --cleanheaders --jobs $jobs ";
	# uncomment below to make separate logs for merge
	# $com .= "2>  $workDir/merge_$caller.log ";
	$com .= "| $bgzip -c -\@12 > $newMerged";
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - starting to merge $caller GVCFs\n";
	system($com) && die "E $0: mergeGvcfs for $caller FAILED: $?";
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - merging $caller GVCFs DONE\n";

	# index:
	# tabix takes a long time on large merged GVCFs but uses a single thread and
	# almost no ressources (RAM or IO): it is cpu-bound but single-threaded
	# -> fork a process for it so we can start working with the next $caller
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - indexing merged $caller GVCF\n";
	$com = "$tabix $newMerged";
	my $pid = fork();
	(defined $pid) ||
	    die "E $0: could not fork for tabix-indexing $caller merged\n";
	if ($pid == 0) {
	    # child: run tabix and exit
	    system($com) && die "E $0: tabix for merged $caller FAILED: $?";
	    # this only killed the child but at least we sent E: message to stderr
	    $now = strftime("%F %T", localtime);
	    warn "I $now: $0 - indexing merged $caller GVCF DONE\n";
	    exit(0);
	}
	else {
	    # parent
	    push(@childrenPids, $pid);
	}
    }
    else {
	# every sample is already in $prevMerged, nothing to do except clean up
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - merging $caller GVCFs not needed, no new samples\n";
	unlink($batchFile) ||
	    die "E $now: $0 - failed to unlink batchFile $batchFile: $!";
    }
}

################################

# wait for tabix processes
foreach my $pid (@childrenPids) {
    waitpid($pid, 0);
}

$now = strftime("%F %T", localtime);
my $mess = "I $now: $0 - ALL DONE!\n";
$mess .= "I $0: Please examine the logs, and if AOK you can remove obsolete merged GVCFs\n";
$mess .= "I $0: ";

my $mirror = &mirror();
($mirror) && ($mess .= "and sync all results to $mirror ");

$mess .= "with the following commands:\n";

$mess .= "cd $dataDir\n";
foreach my $caller (sort(keys %callerDirs)) {
    my $oldestMerged = `ls -rt1 $callerDirs{$caller}->[2]/*.g.vcf.gz | head -n 1`;
    chomp($oldestMerged);
    $mess .= "rm $oldestMerged\n";
    $mess .= "rm $oldestMerged.tbi\n";
}
if ($mirror) {
    $mess .= "rsync -rtvn --delete $bamDir $mirror/$bamDir\n";
    $mess .= "rsync -rtvln --delete $allBamsDir $mirror/$allBamsDir\n";
    $mess .= "rsync -rtvn --delete $gvcfDir $mirror/$gvcfDir\n";
    $mess .= "## redo without -n if AOK:\n";
    $mess .= "rsync -rtv --delete $bamDir $mirror/$bamDir\n";
    $mess .= "rsync -rtvl --delete $allBamsDir $mirror/$allBamsDir\n";
    $mess .= "rsync -rtv --delete $gvcfDir $mirror/$gvcfDir\n";
}

warn "$mess\n";
