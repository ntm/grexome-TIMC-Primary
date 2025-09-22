#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2025
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
# in each job with 16 cores, with 2h time limit.
#
# takes 2 args: $first and $last, the first and last grexomes
# for this run of bam2gvcf_strelka_makeOarJobs.pl
# NOTE: not more than 50*$samplesPerJob at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.


use strict;
use warnings;

# number of grexomes to process in a single OAR job on a 16-core resource in 2h
my $samplesPerJob = 10;

(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50*$samplesPerJob) &&
    die "E: cannot queue more than 50 jobs at a time on OAR/dahu\n";


#################################
# hard-coded stuff

# number of cores/threads (16 on dahu is nice)
my $threads = 16;

# hard-coded dirs:
# log to $logDir
my $logDir = "/home/thierryn/Bam2gvcf_Strelka_stdouterr/";
(-d $logDir) || mkdir($logDir) || 
    die "E: logDir $logDir doesn't exist and can't be mkdir'd\n";
# path to bam2gvcf_strelka.pl
my $bam2gvcf = "/home/thierryn/Bam2gvcf_Strelka/2_bam2gvcf_strelka.pl";
# path/to/latest/configureStrelkaGermlineWorkflow.py
my $strelka = "/home/thierryn/Bam2gvcf_Strelka/strelka-latest/bin/configureStrelkaGermlineWorkflow.py";
# bams are in $inDir
my $inDir = "/bettik/thierryn/BAMs_All_Selected/";
# produce gvcfs in $outDir
my $outDir = "/bettik/thierryn/GVCFs_Strelka_Raw/";
# ref genome
my $genome = "/bettik/thierryn/HumanGenome/hs38DH.fa";
# gzipped and tabix-indexed BED file with chromosomes 1-22, X, Y, M
my $chroms = "/bettik/thierryn/HumanGenome/hs38_chroms.bed.gz";

# oarsub command with params:
# run on my project ngs-timc, ask for 16 cores on 1 node, 2h walltime max
my $oarBase = "oarsub --project ngs-timc -l /nodes=1/core=$threads,walltime=2 ";


#################################
# queue jobs

my @samples = ();
foreach my $gNum ($first..$last) {
    my $grexome = $gNum;
    # left-pad with zeroes to 4 digits
    ($gNum < 10) && ($grexome = "0$grexome");
    ($gNum < 100) && ($grexome = "0$grexome");
    ($gNum < 1000) && ($grexome = "0$grexome");
    $grexome = "grexome$grexome";
    push(@samples, $grexome);
    if ((@samples == $samplesPerJob) || ($gNum == $last)) {
        # choose stdout and stderr filenames
        my $sampsString = $samples[0]."-".$samples[$#samples];
        my $oar = $oarBase."-O $logDir/bam2gvcfStrelka.$sampsString.out -E $logDir/bam2gvcfStrelka.$sampsString.err ";
        $oar .= "\" perl $bam2gvcf --indir $inDir --genome $genome --chroms $chroms --outdir $outDir --samples ";
        $oar .= join(',', @samples);
        $oar .= " --strelka $strelka --jobs $threads --real\"";
        #warn "$oar\n";
        system($oar);
        @samples = ();
    }
}

