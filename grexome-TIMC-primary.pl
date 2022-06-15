#!/usr/bin/perl


# NTM
# 03/02/2021


# This is a wrapper script for the grexome-TIMC primary analysis
# pipeline, starting from "grexomized" FASTQs...
# See $USAGE.

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
# order Raw - Filtered - Merged ($dataDir/$gvcfDir will be prepended soon).
# Here %callerDirs is a template for every variant-caller this pipeline knows,
# --callers will specify the variant-callers to actually use
my %callerDirs = (
    "strelka" => ["GVCFs_Strelka_Raw/","GVCFs_Strelka_Filtered/","GVCFs_Strelka_Filtered_Merged/"],
    "deepvariant" => ["GVCFs_DV_Raw/","GVCFs_DV_Filtered/","GVCFs_DV_Filtered_Merged/"],
    "gatk" => ["GVCFs_GATK_Raw/","GVCFs_GATK_Filtered/","GVCFs_GATK_Filtered_Merged/"],
    "elprep" => ["GVCFs_ElPrep_Raw/","GVCFs_ElPrep_Filtered/","GVCFs_ElPrep_Filtered_Merged/"]);

# merging occurs in batches of up to $mergeBatchSize, and in a second step
# the per-batch merged GVCFs are merged together.
# This is a performance tuning param that doesn't change the final results
my $mergeBatchSize = 25;


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

# comma-separated list of variant-callers to use, default to none ie
# only produce BAMs
my $callers = "";

# workDir must not exist, it will be created and populated
my $workDir;

# type of data among {exome,genome}, used by some variant-callers
my $datatype = 'exome';

# indication of the number of threads / parallel jobs that we can run
my $jobs = 20;

# path+file of the config file holding all install-specific params,
# defaults to the distribution-povided file that you can edit but
# you can also copy it elsewhere and customize it, then use --config
my $config = "$RealBin/grexomeTIMCprim_config.pm";

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "Run the grexome-TIMC primary analysis pipeline, ie start from \"grexomized\" FASTQs and:
- produce BAMs with fastq2bam.pl (trim, align, mark dupes, sort);
- for each specified variant-caller (if any):
     produce individual GVCFs with bam2gvcf_\$caller.pl;
     filter low-quality variant calls with filterBadCalls.pl;
     produce a merged GVCF per variant-caller with mergeGVCFs.pl.

Optionally, if the samplesFile has a \"Sex\" column, also produce qc_sexChroms*.csv 
files in GVCFs_*_Filtered/ , counting HOMO and HET calls on chromosomes X, Y and 16 
(as a control) in every sample, and identifying possible outliers (can indicate 
mis-labeling of samples or other quality issues).

Each step of the pipeline is a stand-alone self-documented script, this is just a wrapper.
For each sample, any step where the result file already exists is skipped.

The \"samples\" must appear in the 'sampleID' column of the samples metadata XLSX file
(provided with --samplesFile).
Default behavior is to process every non-'0' sampleID from this file;
processing can be restricted to specific samples of interest with --SOIs.

The FASTQs must be \"grexomized\", ie for each \$sample we expect a single pair of 
FASTQ files in \$fastqDir (defined in *config.pm), and these files must be named 
\$sample_1.fq.gz and \$sample_2.fq.gz .

BAMs and GVCFs are produced in a hierarchy of subdirs defined at the top of this script,
all within \$dataDir that you should customize (defined in *config.pm). The hierarchy 
can also be modified if desired.
Logs and copies of the metadata are produced in the provided \$workDir .

Every parameter that you could want to customize should be in this script or 
in grexomeTIMCprim_config.pm.

Arguments [defaults] (all can be abbreviated to shortest unambiguous prefixes):
--samplesFile : samples metadata xlsx file, with path
--SOIs : optional, comma-separated list of \"samples of interest\" to process
--callers [default: none, ie only produce BAMs]: comma-separated list of variant-callers to use among ".join(',',keys(%callerDirs))."
--workdir : subdir where logs and workfiles will be created, must not pre-exist
--datatype [$datatype] : type of data among {exome,genome}
--jobs [$jobs] : approximate number of threads/jobs that we should run
--config [defaults to grexomeTIMCprim_config.pm alongside this script] : your customized copy of *config.pm
--help : print this USAGE";

GetOptions ("samplesFile=s" => \$samplesFile,
	    "SOIs=s" => \$SOIs,
	    "callers=s" => \$callers,
	    "workdir=s" => \$workDir,
	    "datatype=s" => \$datatype,
	    "jobs=i" => \$jobs,
	    "config=s" => \$config,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$USAGE\n\n";

($samplesFile) || die "E $0: you must provide a samples metadata file. Try $0 --help\n";
(-f $samplesFile) || die "E $0: the supplied samples metadata file doesn't exist: $samplesFile\n";

{
    # make sure every requested caller is known
    my %requestedCallers;
    foreach my $c (split(/,/,$callers)) {
	# silently allow uppercase or mixed-case
	my $cLow = lc($c);
	defined($callerDirs{$cLow}) ||
	    die "E $0: the provided variant-caller $c is unknown, this pipeline only knows ".
	    join(',',keys(%callerDirs))."\n";
	# OK remember this caller was requested
	$requestedCallers{$cLow} = 1;
    }
    # now delete %callerDirs entries for non-requested callers
    foreach my $c (keys(%callerDirs)) {
	($requestedCallers{$c}) || (delete $callerDirs{$c});
    }
}

# immediately import $config, so we die if file is broken
# if $config doesn't have a path component, prepend ./ to avoid loading the dist version
# (in case the dist was copied into current dir and customized but not renamed)
($config =~ m~/~) || ($config = "./$config");
(-f $config) ||  die "E $0: the supplied config.pm doesn't exist: $config\n";
require($config);
grexomeTIMCprim_config->import( qw(dataDir fastqDir mirror refGenome refGenomeElPrep refGenomeChromsBed 
				   fastTmpPath binPath bwakitPath strelkaBin deepVariantSIF gatkBin elprepBin) );

($workDir) || die "E $0: you must provide a workDir. Try $0 --help\n";
(-e $workDir) && 
    die "E $0: workDir $workDir already exists, remove it or choose another name.\n";

($datatype eq 'exome') || ($datatype eq 'genome') ||
    die "E $0: illegal datatype $datatype, must be among {exome,genome}\n";

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
# parse samples metadata file to grab sampleIDs, limit to --SOIs if specified
# This also serves as an early sanity-check of $samplesFile

# key==existing sample to process
my %samples = ();
#  also test whether we have a "Sex" column in $sampleFile
my $haveSex = 0;

{
    my @parsed = &parseSamples($samplesFile);
    (@parsed == 5) && ($haveSex=1);
    my $s2pathoR = $parsed[0];
    foreach my $s (keys %$s2pathoR) {
	$samples{$s} = 1;
    }
}

if ($SOIs) {
    # make sure every listed sample is in %samples and promote it's value to 2
    foreach my $soi (split(/,/, $SOIs)) {
	($samples{$soi}) ||
	    die "E $0: processing --SOIs: a specified sample $soi does not exist in the samples metadata file\n";
	($samples{$soi} == 1) ||
	    warn "W $0: processing --SOIs: sample $soi was specified twice, is that a typo? Skipping the dupe\n";
	$samples{$soi} = 2;
    }
    # now ignore all other samples
    foreach my $s (keys %samples) {
	if ($samples{$s} != 2) {
	    delete($samples{$s});
	}
    }
}

# exclude any sample that doesn't have FASTQs, but die if called with --SOIs
foreach my $s (sort(keys %samples)) {
    my $f1 = "$fastqDir/${s}_1.fq.gz";
    my $f2 = "$fastqDir/${s}_2.fq.gz";
    if ((! -f $f1) || (! -f $f2)) {
	if ($samples{$s} == 2) {
	    die "E $0: sample $s from --SOIs doesn't have FASTQs (looking for $f1 and $f2)\n";
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
    $com .= "--genome ".&refGenome()." --bwakit ".&bwakitPath()." --threads $jobs --real";
    (&binPath() ne '') && ($com .= " --binpath ".&binPath());
    
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
	my $com = "perl $b2gBin --indir $dataDir/$allBamsDir --samples $samples";
	$com .= " --chroms ".&refGenomeChromsBed();
	$com .= " --outdir $callerWorkDir --jobs $jobs --real";
	#caller-specific args
	if ($caller eq "strelka") {
	    $com .= " --datatype $datatype --genome ".&refGenome();
	    $com .= " --strelka ".&strelkaBin();
	}
	elsif ($caller eq "deepvariant") {
	    $com .= " --datatype $datatype --tmpdir $tmpDir/deepVariant --genome ".&refGenome();
	    $com .= " --deepvariant ".&deepVariantSIF();
	}
	elsif ($caller eq "gatk") {
	    $com .= " --tmpdir $tmpDir/gatk --genome ".&refGenome();
	    $com .= " --gatk ".&gatkBin();
	}
	elsif ($caller eq "elprep") {
	    $com .= " --tmpdir $tmpDir/elprep --logdir $callerWorkDir --genome ".&refGenomeElPrep();
	    # we prefer SFM mode, filter mode is too resource-hungry and not
	    # much faster in our hands. You can test it with "--mode filter"
	    $com .= " --mode sfm";
	    $com .= " --elprep ".&elprepBin();
	}
	else {
	    die "E $0: unknown variant-caller $caller, need to implement caller-specific args";
	}
	
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
	    $com = "perl $RealBin/2_bam2gvcf_strelka_moveGvcfs.pl $callerWorkDir ".$callerDirs{$caller}->[0];
	    system($com) && die "E $0: strelka moveGvcfs FAILED: $?";
	}
	elsif ($caller eq "deepvariant") {
	    # DV log is a single file per sample, not trying to parse it, just moving it with the
	    # GVCF and TBI (and HTML stats file if present) into $gvcfDir subtree and removing
	    # now-empty callerWorkDir:
	    foreach my $s (split(/,/,$samples)) {
		foreach my $file ("$s.g.vcf.gz", "$s.g.vcf.gz.tbi", "$s.log") {
		    move("$callerWorkDir/$file", $callerDirs{$caller}->[0]) ||
			die "E $0: cannot move $callerWorkDir/$file to ".$callerDirs{$caller}->[0]." : $!";
		}
		# hard-coded HTML stats file from DV
		my $htmlStatsDV = "$s.visual_report.html";
		if (-e "$callerWorkDir/$htmlStatsDV") {
		    move("$callerWorkDir/$htmlStatsDV", $callerDirs{$caller}->[0]) ||
			die "E $0: cannot move DV stats $callerWorkDir/$htmlStatsDV to ".$callerDirs{$caller}->[0]." : $!";
		}
	    }
	    rmdir($callerWorkDir) || die "E $0: cannot rmdir $caller callerWorkDir $callerWorkDir: $!";
	}
	elsif ($caller eq "gatk") {
	    # GATK logs are a mess: they seem to adopt a format but then don't respect it,
	    # not even trying to parse it
	    # move GATK GVCFs + TBIs + logs into $gvcfDir subtree and remove now-empty callerWorkDir:
	    foreach my $s (split(/,/,$samples)) {
		foreach my $file ("$s.g.vcf.gz", "$s.g.vcf.gz.tbi", "$s.log") {
		    move("$callerWorkDir/$file", $callerDirs{$caller}->[0]) ||
			die "E $0: cannot move $callerWorkDir/$file to ".$callerDirs{$caller}->[0]." : $!";
		}
	    }
	    rmdir($callerWorkDir) || die "E $0: cannot rmdir gatk callerWorkDir $callerWorkDir: $!";
	}
	elsif ($caller eq "elprep") {
	    # don't yet know what to expect in elprep logs when something goes wrong,
	    # for now just hope that it returns non-zero if there's a problem, we are
	    # keeping the logs so we can investigate if something is fishy...
	    # move elPrep GVCFs into $gvcfDir subtree
	    foreach my $s (split(/,/,$samples)) {
		my $file = "$s.g.vcf.gz";
		move("$callerWorkDir/$file", $callerDirs{$caller}->[0]) ||
		    die "E $0: cannot move $callerWorkDir/$file to ".$callerDirs{$caller}->[0]." : $!";
	    }
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
	# samples to filter: those without a filtered GVCF
	my $gvcf = $callerDirs{$caller}->[1]."/$s.filtered.g.vcf.gz";
	if (! -e $gvcf) {
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
	if (! -e "$gvcf.tbi") {
	    # index filtered GVCF
	    my $com = "$tabix -p vcf $gvcf";
	    system($com) && die "E $0: tabix-index for individual GVCF $gvcf FAILED: $?";
	}
    }
    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - filtering $caller GVCFs DONE\n";

    ################################
    # create/update QC file that counts HET/HV calls on sex chromosomes,
    # this requires knowing if samples are M or F
    if ($haveSex) {
	# name of file to create/update
	my $qcFile = $callerDirs{$caller}->[1]."/qc_sexChroms_$caller.csv";

	my $com = "perl $RealBin/0_qc_checkSexChroms.pl  --samplesFile=$samplesFile";
	$com .= " --indir=".$callerDirs{$caller}->[1];
	$com .= " --tabix=$tabix";
	# backup and use previous version (if any)
	if (-e $qcFile) {
	    my $qcPrev = $qcFile;
	    ($qcPrev =~ s/\.csv$/_prev.csv/) ||
		die "E $0: cannot subst csv in qcFile $qcFile, WTF?!\n";
	    move($qcFile,$qcPrev) ||
		die "E $0: qc file $qcFile exists but can't be moved to $qcPrev: $!";
	    # we have a backup => can --force
	    $com .= " --prevQC=$qcPrev --force";
	}
	$com .= " > $qcFile";
	system($com) && die "E $0: qc_checkSexChroms for $caller FAILED: $?";
    }
    else {
	$now = strftime("%F %T", localtime);
	warn "W $now: $0 - no Sex column in metadata, skipping qc_checkSexChroms step\n";
    }

    ################################
    # merge new GVCFs with the most recent previous merged if one existed

    # samples in prevMerged, to avoid dupes
    my %samplesPrev;

    # array of filenames (with path) of the new filtered GVCFs
    my @newToMerge = ();
    
    # we want to merge the new GVCFs with the most recent previous merged,
    # if there was one. code is a bit ugly but functional
    my $prevMerged = `ls -rt1 $callerDirs{$caller}->[2]/*.g.vcf.gz 2> /dev/null | tail -n 1`;
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
    }
    
    foreach my $s (sort(keys %samples)) {
	# only merge $s if it's not already in prevMerged
	if (! $samplesPrev{$s}) {
	    push(@newToMerge, $callerDirs{$caller}->[1]."/$s.filtered.g.vcf.gz");
	}
    }

    # if we have more than $mergeBatchSize new GVCFs, merge them in batches
    if (@newToMerge >= $mergeBatchSize) {
	my $batchNum = 0;
	my @batchFiles;
	my $batchFH;
	foreach my $i (0..$#newToMerge) {
	    if ($i % $mergeBatchSize == 0) {
		$batchNum++;
		($batchFH) && close($batchFH);
		my $batchFile = "$workDir/batchFile_${caller}_BATCH$batchNum.txt";
		push(@batchFiles, $batchFile);
		open($batchFH, ">$batchFile") ||
		    die "E $0: cannot create $caller batchFile $batchFile: $!\n";
	    }
	    print $batchFH $newToMerge[$i]."\n";
	}
	# close last batchFile
	close($batchFH);

	# merge each batch
	my @batchGVCFs;
	foreach my $b (1..$batchNum) {
	    my $newMerged = $callerDirs{$caller}->[2]."/grexomes_${caller}_merged_${date}_BATCH$b.g.vcf.gz";
	    (-e $newMerged) &&
		die "E $0: want to batchwise-merge GVCFs but newMerged already exists: $newMerged\n";
	    push(@batchGVCFs, $newMerged);
	    my $batchFile = $batchFiles[$b-1];
	    # -> merge:
	    my $com = "perl $RealBin/4_mergeGVCFs.pl --filelist $batchFile --tmpdir $tmpDir/Merge --cleanheaders --jobs $jobs ";
	    # uncomment below to make separate logs for merge
	    # $com .= " 2> $workDir/merge_${caller}_BATCH${batch}.log ";
	    # NOTE we are piping to bgzip -@8 , so we exceed $jobs quite a bit.
	    # if this turns into a problem we can tune it down (eg $jobsFilter)
	    $com .= "| $bgzip -c -\@8 > $newMerged";
	    $now = strftime("%F %T", localtime);
	    warn "I $now: $0 - starting to batchwise-merge $caller GVCFs batch $b\n";
	    system($com) && die "E $0: batchwise-mergeGvcfs for $caller batch $b FAILED: $?";
	    $now = strftime("%F %T", localtime);
	    warn "I $now: $0 - batchwise-merging $caller GVCFs batch $b DONE\n";
 	}

	# now replace individual files in @newToMerge by the new batchwise GVCFs,
	# so the final merge works the same whether we did batchwise merging or not
	@newToMerge = @batchGVCFs;
    }

    # only merge if there's at least one new sample
    if (@newToMerge) {
	my $newMerged = $callerDirs{$caller}->[2]."/grexomes_${caller}_merged_$date.g.vcf.gz";
	(-e $newMerged) &&
	    die "E $0: want to merge GVCFs but newMerged already exists: $newMerged\n";
	# make batchfile with list of GVCFs to merge
	my $batchFile = "$workDir/batchFile_$caller.txt";
	open(my $bf, ">$batchFile") ||
	    die "E $0: cannot create $caller batchFile $batchFile: $!\n";
	($prevMerged) && (print $bf "$prevMerged\n");
	foreach my $new (@newToMerge) {
	    print $bf "$new\n";
	}
	close($bf);

	# -> merge:
	# NOTE we are piping to bgzip -@12 , so we exceed $jobs quite a bit.
	# if this turns into a problem we can tune it down (eg $jobsFilter)
	my $com = "perl $RealBin/4_mergeGVCFs.pl --filelist $batchFile --tmpdir $tmpDir/Merge --cleanheaders --jobs $jobs ";
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
	$com = "$tabix -p vcf $newMerged";
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
	# every sample is already in $prevMerged, nothing to do
	$now = strftime("%F %T", localtime);
	warn "I $now: $0 - merging $caller GVCFs not needed, no new samples\n";
    }
}

################################

# wait for tabix processes
foreach my $pid (@childrenPids) {
    waitpid($pid, 0);
}

$now = strftime("%F %T", localtime);
my $mess = "I $now: $0 - ALL DONE!\n";
$mess .= "I $0: Please examine the logs.\n";

$mess .= "I $0: If AOK you can remove obsolete merged GVCFs if any (keeping 2 for each caller)\n";
$mess .= "I $0: ";

my $mirror = &mirror();
($mirror) && ($mess .= "and sync all results to $mirror ");

$mess .= "with the following commands:\n";

$mess .= "cd $dataDir\n";
foreach my $caller (sort(keys %callerDirs)) {
    my @GvcfsToRm = split("\n",`ls -rt1 $callerDirs{$caller}->[2]/*.g.vcf.gz | grep -v _BATCH`);
    # keep two most recent
    pop(@GvcfsToRm);
    pop(@GvcfsToRm);
    my @TBIsToRm;
    foreach my $gvcf (@GvcfsToRm) {
	(-e "$gvcf.tbi") && push(@TBIsToRm,"$gvcf.tbi");
    }
    # batchwise merged GVCFs were temp, remove them all
    push(@GvcfsToRm, split("\n",`ls -rt1 $callerDirs{$caller}->[2]/*.g.vcf.gz | grep _BATCH`));
    (@GvcfsToRm) && ($mess .= "rm ".join(" ", @GvcfsToRm)."\n");
    (@TBIsToRm) && ($mess .= "rm ".join(" ", @TBIsToRm)."\n");
}
if ($mirror) {
    $mess .= "rsync -rtvn --delete $bamDir $mirror/$bamDir\n";
    $mess .= "rsync -rtvln --delete $allBamsDir $mirror/$allBamsDir\n";
    $mess .= "rsync -rtvln --delete $gvcfDir $mirror/$gvcfDir\n";
    $mess .= "## redo without -n if AOK:\n";
    $mess .= "rsync -rtv --delete $bamDir $mirror/$bamDir\n";
    $mess .= "rsync -rtvl --delete $allBamsDir $mirror/$allBamsDir\n";
    $mess .= "rsync -rtvl --delete $gvcfDir $mirror/$gvcfDir\n";
}

warn "$mess\n";
