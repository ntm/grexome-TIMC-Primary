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


# 19/06/2019
# NTM
#
# Submit fastq2bam jobs on a cluster with OAR.
# Not sure if this script can be useful outside our organization:
# it would need adapting for a different batch scheduler and/or
# cluster, paths and filename patterns are hard-coded (eg samples
# are called grexome\d\d\d\d), etc... 
#
# takes 2 args: $first and $last, these are ints and assumes that
# all samples are named grexome\d\d\d\d , the 4-digit integer will 
# range from $first to $last for this run of fastq2bam_makeOarJobs.pl.
# NOTE: not more than 50 jobs at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.

use strict;
use warnings;


# number of samples to process in a single OAR job
my $samplesPerJob = 40;

(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50 * $samplesPerJob) &&
    die "E: cannot queue more than 50 jobs == ".50 * $samplesPerJob." samples ($samplesPerJob per job) at a time on OAR/dahu\n";

#################################
# hard-coded stuff

# hard-coded dirs:
# log to $logDir
my $logDir = "/home/thierryn/Fastq2Bam_stdouterr/";
(-d $logDir) || mkdir($logDir) || 
    die "E: logDir $logDir doesn't exist and can't be mkdir'd\n";
# binaries are in $binDir
my $binDir = "/home/thierryn/Software/Fastq2Bam/";
# fastqs are in $inDir
my $inDir = "/bettik/thierryn/FASTQs_All_Grexomized/";
# produce bams in $outDir
my $outDir = "/bettik/thierryn/BAMs_grexome/";
# ref genome
my $genome = "/bettik/thierryn/HumanGenome/hs38DH.fa";

# OAR resources:
# number of cores/threads (16 on dahu was nice, let's try 48 (ie two "numa nodes") on kraken),
# $threads and $resources must correspond, this depends on the architecture of the cluster nodes
my $threads = 48;
my $resources = "/nodes=1/cpu=1/numa=2";
# alternative: $resources = "/core=$threads" and then rely on the [ANTIFRAG] rule
# to run the job on a single CPU and full numa nodes

# oarsub command with params:
# run on my project ngs-timc, ask for $resources and 7h walltime max (must be
# enough for $samplesPerJob samples, 6h was barely enough in 21/01/26 run, 
# 7h should be fine for any future jobs on the same cluster)
my $oarBase = "oarsub --project ngs-timc -l $resources,walltime=7 ";


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
        my $oar = $oarBase."-O $logDir/fastq2bam.$sampsString.$threads.out -E $logDir/fastq2bam.$sampsString.$threads.err ";
        $oar .= "\"perl $binDir/1_fastq2bam.pl --out $outDir --in $inDir --tmp /var/tmp/$sampsString --samples ";
        $oar .= join(',', @samples);
        # no more "--bwakit $binDir" , it is super slow and doesn't seem useful for us
        $oar .= " --bin $binDir --threads $threads --genome $genome --real\"";
        system($oar);
        @samples = ();
    }
}

