#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2022
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


# 08/06/2022
# NTM

# Call variants on BAMs and produce GVCF files using deepVariant.
# Currently requires singularity/apptainer, though the code could be easily
# adapted to run a natively compiled version (or docker image).
#
# See $USAGE for arguments.

use strict;
use warnings;
use Getopt::Long;
use POSIX qw(strftime);
use File::Copy qw(move);
use File::Basename qw(basename fileparse);
use Cwd qw(abs_path);
use File::Path qw(remove_tree);


# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);

#############################################
## options / params from the command-line

# subdir where BAMs can be found
my $inDir;

# comma-separated list of samples (FASTQs) to process (required)
my $samples = '';

# path+filename of ref genome, currently we recommend the full GRCh38 with
# decoy+alts+unmapped, as produced by Heng Li's run-gen-ref (from bwa-kit)
my $refGenome;

# bgzipped and tabix-indexed BED defining regions where variants should be called,
# any other genomic region is ignored
my $chromsBed;

# dir where GVCFs will be created
my $outDir;

# tmpDir, must not pre-exist and will be rm'd, faster is better 
# (ramdisk or at least SSD)
my $tmpDir;

# path+name of deepvariant singularity image, eg "deepvariant_1.4.0.sif"
my $deepVariant = "";

# type of data: exome or genome
my $datatype = "exome";

# number of cores to use (only for make_examples step of DV)
my $jobs = 16;

# $real: if not true don't actually process anything, just print INFO 
# messages on what would be done. Default to false
my $real = '';

# help: if true just print $USAGE and exit
my $help = '';

my $USAGE = "
Arguments (all can be abbreviated to shortest unambiguous prefixes):
--indir : subdir containing the BAMs
--samples : comma-separated list of sampleIDs to process, for each sample we expect
	  [sample].bam and [sample].bam.bai files in indir
--genome : ref genome fasta, with path
--chroms : optional, if provided it must be a bgzipped and tabix-indexed BED file
	   defining regions where variants should be called
--outdir : dir where GVCF files will be created
--tmpdir : subdir where tmp files will be created, must not pre-exist and will be removed after execution
--deepvariant : path+name of deepvariant singularity image (SIF format)
--datatype [$datatype] : type of data, among {exome,genome}
--jobs [$jobs] : number of cores that DV can use
--real : actually do the work, otherwise this is a dry run
--help : print this USAGE";


GetOptions ("indir=s" => \$inDir,
	    "samples=s" => \$samples,
	    "genome=s" => \$refGenome,
	    "chroms=s" => \$chromsBed,
	    "outdir=s" => \$outDir,
	    "tmpdir=s" => \$tmpDir,
	    "deepvariant=s" => \$deepVariant,
	    "datatype=s" => \$datatype,
	    "jobs=i" => \$jobs, 
	    "real" => \$real,
	    "help" => \$help)
    or die("E $0: Error in command line arguments\n$USAGE\n");

# make sure required options were provided and sanity check them
($help) && die "$0 $USAGE\n\n";

($inDir) ||
    die "E $0: you MUST provide --indir where BAMs can be found\n$USAGE\n";
(-d $inDir) ||
    die "E $0: inDir specified is not a folder!";

# save samples in %samples to detect duplicates and allow sorting
my %samples;
foreach my $sample (split(/,/, $samples)) {
    if ($samples{$sample}) {
	warn "W $0: sample $sample was specified twice, is that a typo? Ignoring the dupe\n";
	next;
    }
    $samples{$sample} = 1;
}

($refGenome) || die "E $0: you must provide a ref genome fasta file\n";
(-f $refGenome) || die "E $0: provided genome fasta file doesn't exist\n";

if ($chromsBed) {
    (-f $chromsBed) || die "E $0: provided --chroms file doesn't exist\n";
    (-f "$chromsBed.tbi") || (-f "$chromsBed.csi") ||
	die "E $0: can't find tabix index for provided --chroms file\n";
}

# make sure DV image exists and can be run with singularity
($deepVariant) || die "E $0: you must provide a deepVariant singularity image file\n";
(-f $deepVariant) || die "E $0: provided deepVariant SIF file doesn't exist\n";
$deepVariant .= " /opt/deepvariant/bin/run_deepvariant";
(`singularity run $deepVariant --version` =~ /^DeepVariant version/) ||
    die "E $0: provided DV image broken / mis-behaving, can't find DV version with:\n".
    "\tsingularity run $deepVariant --version\n";

($datatype eq 'exome') || ($datatype eq 'genome') ||
    die "E $0: illegal datatype $datatype, must be among {exome,genome}\n";

($jobs > 0) ||
    die "E $0: called with jobs=$jobs but we need at least one thread\n";

($outDir) || 
    die "E $0: you MUST specify --outdir where GVCFs will be created\n$USAGE\n";
(-d $outDir) || (mkdir($outDir)) || 
    die "E $0: outDir $outDir doesn't exist as a dir and can't be created\n";

($tmpDir) || die "E $0: you must provide a tmpDir\n";
(-e $tmpDir) && 
    die "E $0: called with tmpdir=$tmpDir but it already exists, remove it or choose another name.\n";
mkdir($tmpDir) || die "E $0: cannot mkdir tmpDir $tmpDir\n";


#############################################
## process each sample of interest

my $now = strftime("%F %T", localtime);
warn "I $now: $0 - STARTING TO WORK\n";

if ($chromsBed) {
    # DV needs a gunzipped $chromsBed, make a copy in $tmpDir
    if (system("gunzip -c $chromsBed > $tmpDir/chroms.bed")) {
	# non-zero status, clean up and die
	remove_tree($tmpDir);
	die "E $0: need gunzipped chromsBed but failure with: gunzip -c $chromsBed > $tmpDir/chroms.bed\n";
    }
}

foreach my $sample (sort keys(%samples)) {
    # make sure we have bam and bai files for $sample, otherwise skip
    my $bam = "$inDir/$sample.bam";
    ((-e $bam) && (-e "$bam.bai")) || 
	((warn "W $0: no BAM or BAI for $sample in inDir $inDir, skipping $sample\n") && next);

    #############################################
    ## build the DV command-line

    # need to bind the symlink-resolved dir containing $bam, as well as $outdir,
    # $tmpdir, and the symlink-resolved dirs containing $refGenome, to absolute
    # paths in the singularity container. We will bind to:
    my ($singIn, $singOut, $singTmp, $singRefGen, $singChrBed) =
	("/IN/", "/OUT/", "/TMP/", "/REFGEN/", "/CHRBED/");

    my $bamResolved = abs_path($bam);
    my ($bamResFile,$bamResDir) = fileparse($bamResolved);
    my $bindings = "$bamResDir:$singIn";
    $bindings .= ",$outDir:$singOut,$tmpDir:$singTmp";
    my $refGenResolved = abs_path($refGenome);
    my ($refGenResFile,$refGenResDir) = fileparse($refGenResolved);
    $bindings .= ",$refGenResDir:$singRefGen";
    my $chrResolved = abs_path($chromsBed);

    # override env $TMPDIR (DV writes stuff there) and silence perl warnings due
    # to missing locales in the DV docker/singularity image, but do it all in a 'bash -c " ('
    # so we can capture stderr in OAR (remeber to close the ') "' after redirecting stderr)
    my $com = "bash -c \" ( TMPDIR=$singTmp LC_ALL=C ";
    $com .= " singularity run --bind $bindings $deepVariant";

    $com .= " --reads=$singIn/$bamResFile";
    # gvcf to produce
    my $gvcf = "${sample}.g.vcf.gz";
    (-e "$outDir/$gvcf") && 
	(warn "W $0: GVCF for $sample already exists in outDir $outDir, skipping $sample\n") && next;
    $com .= " --output_gvcf=$singOut/$gvcf";

    # vcf to produce: we only want GVCFs but DV insists on also creating VCFs,
    # we'll remove them when done
    my $vcf = "$singTmp/$sample.vcf.gz";
    $com .= " --output_vcf=$vcf";

    if ($datatype eq 'exome') { $com .= " --model_type=WES"; }
    else { $com .= " --model_type=WGS"; }
    $com .= " --ref=$singRefGen/$refGenResFile";
    ($chromsBed) && ($com .= " --regions=$singTmp/chroms.bed");
    $com .= " --intermediate_results_dir=$singTmp/intermediate/";
    $com .= " --num_shards=$jobs";
    # keep HTML stats file for now, to disable uncomment:
    # $com .= " --novcf_stats_report"; # disable HTML stats logfile
    ($real) || ($com .= " --dry_run=true");
    # unused: --logging_dir

    
    # DV logging: one file per sample
    my $log = "${sample}.log";
    # need to close parenthesis and quote (for capturing stderr)
    $com .= " &> $outDir/$log ) \" ";

    $now = strftime("%F %T", localtime);
    warn "I $now: $0 - running deepVariant for $sample\n";
    if (system($com)) {
	# non-zero status, clean up and die
	remove_tree($tmpDir);
	$now = strftime("%F %T", localtime);
	(-e $tmpDir) && 
	    warn "E $now: $0 - running deepVariant FAILED but cannot rmdir tmpDir $tmpDir, why?\n";
        die "E $now: $0 - running deepVariant for $sample FAILED ($?)! INSPECT THE LOGFILE $outDir/$log\n";
    }
    else {
	# DV completed successfully.
	# If we didn't disable HTML stats, the filename can't be specified
	# so we need to rename it or it'll get squashed by the next sample
	my $htmlStatsDV = "$outDir/output.visual_report.html";
	if (-e "$htmlStatsDV") {
	    my $htmlStatsNew = "$outDir/$sample.visual_report.html";
	    (-e "$htmlStatsNew") &&
		warn "W $0: DV for $sample successful but stats file $htmlStatsNew pre-exists, I will clobber it\n";
	    move("$htmlStatsDV", "$htmlStatsNew") ||
		die "E $0: move failed for $htmlStatsDV to $htmlStatsNew : $!";
	}
    }
}

$now = strftime("%F %T", localtime);
remove_tree($tmpDir);
(-e $tmpDir) && 
    warn "E $now: $0 - all done but cannot rmdir tmpDir $tmpDir, why?\n";
warn "I $now: $0 - ALL DONE\n";
