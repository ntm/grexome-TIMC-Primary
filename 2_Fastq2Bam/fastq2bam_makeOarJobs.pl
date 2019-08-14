#!/usr/bin/perl

# 19/06/2019
# NTM
#
# submit fastq2bam jobs on dahu with oar.
# I don't want to modify fastq2bam.pl (eg it's good that all the
# sanity checks are performed on each node / by each job), so I
# am just calling fastq2bam.pl many times on individual samples.


use strict;
use warnings;

# takes 2 args: $first and $last, the first and last grexomes
# for this run of fastq2bam_makeOarJobs.pl
# NOTE: not more than 50 at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.
(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50) &&
    die "E: cannot queue more than 50 jobs at a time on OAR/dahu\n";

#################################
# hard-coded stuff

# number of cores/threads (16 on dahu is nice)
my $threads = 16;

# hard-coded dirs:
# log to $logDir
my $logDir = "~/Fastq2Bam_Dahu_stdouterr/";
(-d $logDir) || mkdir($logDir) || 
    die "E: logDir $logDir doesn't exist and can't be mkdir'd\n";
# binaries are in $binDir
my $binDir = "~/Fastq2Bam_PackagedWithBinaries_Centos7/";
# fastqs are in $inDir
my $inDir = "/bettik/nthierry/FASTQs_All_Grexomized/";
# produce bams in $outDir
my $outDir = "/bettik/nthierry/BAMs_grexome_NTM/BAMs_NTM_Dahu/";

# oarsub command with params:
# run on my project ngs-timc, ask for 16 cores on 1 node, 2h walltime max
my $oarBase = "oarsub --project ngs-timc -l /nodes=1/core=$threads,walltime=2 ";


#################################
# queue jobs


foreach my $gNum ($first..$last) {
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";
    # choose stdout and stderr filenames
    my $oar = $oarBase."-O $logDir/fastq2bam.$grexome.out -E $logDir/fastq2bam.$grexome.err ";
    $oar .= "\"perl $binDir/fastq2bam.pl --out $outDir --in $inDir --first $gNum --last $gNum --bin $binDir --bwakit $binDir --threads $threads --real\"";

    system($oar);

}

