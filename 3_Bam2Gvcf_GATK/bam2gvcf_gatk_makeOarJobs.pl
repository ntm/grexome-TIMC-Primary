#!/usr/bin/perl


# 08/10/2020
# NTM
#
# submit bam2gvcf_gatk.pl jobs on dahu with oar.
# Not sure if this script can be useful outside our organization:
# it would need adapting for a different batch scheduler and/or
# cluster, paths and filename patterns are hard-coded (eg samples
# are called grexome\d\d\d\d), etc... 
#
# GATK is executed via singularity (calling gatk directly results
# in death with a cryptic error message).
#
# Takes 2 args: $first and $last, the first and last grexomes (samples)
# for this run of bam2gvcf_gatk_makeOarJobs.pl
# NOTE: not more than 50*$grexomesPerJob at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.


use strict;
use warnings;

# max number of grexomes to process in a single OAR job of 24h, 
# assume 8h for 24 samples in parallel -> up to 72 should be OK
my $grexomesPerJob = 72;

(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50*$grexomesPerJob) &&
    die "E: cannot queue more than 50 jobs at a time on OAR/dahu\n";


#################################
# hard-coded stuff

# number of parallel jobs on one node
my $jobs = 24;

# hard-coded dirs:
# log to $logDir
my $logDir = "/bettik/nthierry/ProcessBams_GATK_2010_dahu_logs/";
# path to bam2gvcf_gatk.pl
my $bam2gvcf = "~/Bam2gvcf_GATK_PackagedWithBinaries_Centos7/bam2gvcf_gatk.pl";
# path/to/config.pm
my $config = "~/Bam2gvcf_GATK_PackagedWithBinaries_Centos7/grexomeTIMCprim_config.pm";
# bams are in $inDir
my $inDir = "/bettik/nthierry/BAMs_All_Selected/";
# produce gvcfs in $outDir
my $outDir = "/bettik/nthierry/ProcessBams_GATK_2010_dahu/";

# path/to/latest/gatk : doesn't work, gatk dies with cryptic message
#my $gatk = "~/Bam2gvcf_GATK_PackagedWithBinaries_Centos7/gatk-latest/gatk";
# -> instead we use a singularity image
my $gatk = 'singularity exec';
# bind /bettik/nthierry/ (ie make it rw-accessible from within the container),
# /home/nthierry/ is bound by default
$gatk .= ' --bind /bettik/nthierry/';
# path/to/image
##$gatk .= ' ~/Bam2gvcf_GATK_PackagedWithBinaries_Centos7/gatk-latest.sif';
# running from the sif doesn't work, trying a sandbox
$gatk .= ' ~/gatk_4.1.8.1_SANDBOX';

# run gatk in bash -c with leading "( so we can capture stderr... $bam2gvcf must
# close the paren and quote after redirecting stderr
$gatk .= ' bash -c \" ( gatk';

# oarsub command with params: run on my project ngs-timc, and...
my $oarBase = "oarsub --project ngs-timc";
## INITIAL BIG JOB TO PROCESS GREXOMES 50-489: ask for 1 full node, 24h walltime
##$oarBase .= " -l /nodes=1,walltime=24 ";
# for a single sample: 4 cores on 1 node for 12h:
$oarBase .= " -l /nodes=1/core=4,walltime=12 ";

my $grex1 = $first;
my $grex2 = $first + $grexomesPerJob - 1;
while ($grex1 <= $last) {
    ($grex2 <= $last) || ($grex2 = $last);
    # build comma-separated list of samples
    my $samples = $grex1;
    # left-pad with zeroes to 4 digits
    ($samples < 10) && ($samples = "0$samples");
    ($samples < 100) && ($samples = "0$samples");
    ($samples < 1000) && ($samples = "0$samples");
    $samples = "grexome$samples";
    foreach my $s ($grex1+1..$grex2) {
	($s < 10) && ($s = "0$s");
	($s < 100) && ($s = "0$s");
	($s < 1000) && ($s = "0$s");
	$samples .= ",grexome$s";
    }
    my $oar = $oarBase." -O $logDir/bam2gvcf.$grex1-$grex2.out -E $logDir/bam2gvcf.$grex1-$grex2.err ";
    $oar .= "\" perl $bam2gvcf --indir $inDir --samples $samples --outdir $outDir --gatk \'$gatk\' --config $config --jobs $jobs --real\"";

    #warn "$oar\n";
    system($oar);

    $grex1 = $grex2 + 1;
    $grex2 = $grex1 + $grexomesPerJob - 1;
}
