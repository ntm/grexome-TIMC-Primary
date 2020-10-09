#!/usr/bin/perl


# 08/10/2020
# NTM
#
# submit bam2gvcf_gatk.pl jobs on dahu with oar.
# Given timings on luxor, I should be fine processing 10 grexomes
# in each job, with 2h time limit.
#
# takes 2 args: $first and $last, the first and last grexomes
# for this run of bam2gvcf_gatk_makeOarJobs.pl
# NOTE: not more than 50*$grexomesPerJob at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.


use strict;
use warnings;

# number of grexomes to process in a single OAR job of 24h, 
# assume 6h for 24 samples in parallel -> up to 96 could be OK
my $grexomesPerJob = 82;

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
# path/to/latest/gatk
my $gatk = "~/Bam2gvcf_GATK_PackagedWithBinaries_Centos7/gatk-latest/gatk";
# path/to/config.pm
my $config = "~/Bam2gvcf_GATK_PackagedWithBinaries_Centos7/grexomeTIMCprim_config.pm";
# bams are in $inDir
my $inDir = "/bettik/nthierry/BAMs_All_Selected/";
# produce gvcfs in $outDir
my $outDir = "/bettik/nthierry/ProcessBams_GATK_2010_dahu/";

# oarsub command with params:
# run on my project ngs-timc, ask for 1 full node, 24h walltime max
my $oarBase = "oarsub --project ngs-timc -l /nodes=1,walltime=24 ";

my $grex1 = $first;
my $grex2 = $first + $grexomesPerJob - 1;
while ($grex1 <= $last) {
    ($grex2 <= $last) || ($grex2 = $last);
    my $oar = $oarBase."-O $logDir/bam2gvcf.$grex1-$grex2.out -E $logDir/bam2gvcf.$grex1-$grex2.err ";
    $oar .= "\"perl $bam2gvcf --indir $inDir --outdir $outDir --gatk $gatk --first $grex1 --last $grex2 --config $config --jobs $jobs --real\"";

    system($oar);
    
    $grex1 = $grex2 + 1;
    $grex2 = $grex1 + $grexomesPerJob - 1;
}

