#!/usr/bin/perl


# 22/06/2019
# NTM
#
# submit bam2gvcf_strelka.pl jobs on dahu with oar.
# Not sure if this script can be useful outside our organization:
# it would need adapting for a different batch scheduler and/or
# cluster, paths and filename patterns are hard-coded (eg samples
# are called grexome\d\d\d\d), etc... 
#
# Given timings on luxor, I should be fine processing 10 grexomes
# in each job, with 2h time limit.
#
# takes 2 args: $first and $last, the first and last grexomes
# for this run of bam2gvcf_strelka_makeOarJobs.pl
# NOTE: not more than 50*$grexomesPerJob at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.


use strict;
use warnings;

# number of grexomes to process in a single OAR job of 2h, 10 should be fine
my $grexomesPerJob = 10;

(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50*$grexomesPerJob) &&
    die "E: cannot queue more than 50 jobs at a time on OAR/dahu\n";


#################################
# hard-coded stuff

# number of cores/threads (16 on dahu is nice)
my $threads = 16;

# hard-coded dirs:
# log to $logDir
my $logDir = "/bettik/thierryn/ProcessBams_Strelka_1906_dahu_logs/";
# path to bam2gvcf_strelka.pl
my $bam2gvcf = "~/Bam2gvcf_Strelka_PackagedWithBinaries_Centos7/2_bam2gvcf_strelka.pl";
# path/to/latest/configureStrelkaGermlineWorkflow.py
my $strelka = "~/Bam2gvcf_Strelka_PackagedWithBinaries_Centos7/strelka-latest/bin/configureStrelkaGermlineWorkflow.py";
# bams are in $inDir
my $inDir = "/bettik/thierryn/BAMs_All_Selected/";
# produce gvcfs in $outDir
my $outDir = "/bettik/thierryn/ProcessBams_Strelka_1906_dahu/";

# oarsub command with params:
# run on my project ngs-timc, ask for 16 cores on 1 node, 2h walltime max
my $oarBase = "oarsub --project ngs-timc -l /nodes=1/core=$threads,walltime=2 ";

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
    my $oar = $oarBase."-O $logDir/bam2gvcf.$grex1-$grex2.out -E $logDir/bam2gvcf.$grex1-$grex2.err ";
    $oar .= "\"perl $bam2gvcf --in $inDir --samples $samples --out $outDir --strelka $strelka --jobs $threads --real\"";
    system($oar);
    
    $grex1 = $grex2 + 1;
    $grex2 = $grex1 + $grexomesPerJob - 1;
}

