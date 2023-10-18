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


# 19/06/2019
# NTM
#
# submit fastq2bam jobs on the dahu cluster with OAR.
# Not sure if this script can be useful outside our organization:
# it would need adapting for a different batch scheduler and/or
# cluster, paths and filename patterns are hard-coded (eg samples
# are called grexome\d\d\d\d), etc... 

use strict;
use warnings;

# takes 2 args: $first and $last, these are ints and assumes that
# all samples are named grexome\d\d\d\d , the 4-digit integer will 
# range from $first to $last for this run of fastq2bam_makeOarJobs.pl.
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
my $logDir = "/home/thierryn/Fastq2Bam_Dahu_stdouterr/";
(-d $logDir) || mkdir($logDir) || 
    die "E: logDir $logDir doesn't exist and can't be mkdir'd\n";
# binaries are in $binDir
my $binDir = "/home/thierryn/Fastq2Bam_PackagedWithBinaries/";
# fastqs are in $inDir
my $inDir = "/bettik/thierryn/FASTQs_All_Grexomized/";
# produce bams in $outDir
my $outDir = "/bettik/thierryn/BAMs_grexome_NTM/BAMs_NTM_Dahu/";
# ref genome
my $genome = "/bettik/thierryn/HumanGenome/hs38DH.fa";

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
    $oar .= "\"perl $binDir/1_fastq2bam.pl --out $outDir --in $inDir --samples $grexome --bin $binDir --bwakit $binDir --threads $threads --genome $genome -real\"";

    system($oar);

}

