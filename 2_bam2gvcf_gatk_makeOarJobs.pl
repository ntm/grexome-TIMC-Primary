#!/usr/bin/perl


############################################################################################
# Copyright (C) Nicolas Thierry-Mieg, 2019-2024
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


# 08/10/2020
# NTM
#
# submit bam2gvcf_gatk.pl jobs on the dahu cluster with OAR.
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
# NOTE: not more than 50*$samplesPerJob at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.

use strict;
use warnings;


# max number of grexomes to process in a single OAR job:
# trying to run 8 samples in parallel on a 16-core resource in 8h
my $samplesPerJob = 8;

(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50*$samplesPerJob) &&
    die "E: cannot queue more than 50 jobs == ".50 * $samplesPerJob." samples ($samplesPerJob per job) at a time on OAR/dahu\n";


#################################
# hard-coded stuff

# number of cores/threads (16 on dahu is nice)
my $jobs = 16;

# hard-coded dirs:
# log to $logDir
my $logDir = "/home/thierryn/Bam2gvcf_GATK_stdouterr/";
(-d $logDir) || mkdir($logDir) || 
    die "E: logDir $logDir doesn't exist and can't be mkdir'd\n";
# path to bam2gvcf_gatk.pl
my $bam2gvcf = "/home/thierryn/Bam2gvcf_GATK/2_bam2gvcf_gatk.pl";
# bams are in $inDir
my $inDir = "/bettik/thierryn/BAMs_All_Selected/";
# produce gvcfs in $outDir
my $outDir = "/bettik/thierryn/GVCFs_GATK_Raw/";
# ref genome
my $genome = "/bettik/thierryn/HumanGenome/hs38DH.fa";
# gzipped and tabix-indexed BED file with chromosomes 1-22, X, Y, M
my $chroms = "/bettik/thierryn/HumanGenome/hs38_chroms.bed.gz";

# path/to/latest/gatk : doesn't work, gatk dies with cryptic message
#my $gatk = "/home/nthierry/Bam2gvcf_GATK/gatk-latest/gatk";
# -> instead we use a singularity image
my $gatk = 'singularity exec';
# bind /bettik/thierryn/ (ie make it rw-accessible from within the container),
# /home/thierryn/ is bound by default
$gatk .= ' --bind /bettik/thierryn/';
# path/to/image
$gatk .= ' /home/thierryn/Bam2gvcf_GATK/gatk-latest.sif';
# running from the sif doesn't work, trying a sandbox
# $gatk .= ' ~/gatk_4.1.8.1_SANDBOX';

# run gatk in bash -c with leading "( so we can capture stderr... $bam2gvcf must
# close the paren and quote after redirecting stderr
$gatk .= ' bash -c \" ( gatk';

# oarsub command with params: run on my project ngs-timc, and...
my $oarBase = "oarsub --project ngs-timc";
# many samples: ask for 16 cores on 1 node, 8h walltime (should be enough for 8 samples)
$oarBase .= " -l /nodes=1/core=$jobs,walltime=8 ";
# for a single sample: 4 cores on 1 node for 12h:
## $oarBase .= " -l /nodes=1/core=4,walltime=12 ";


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
	my $oar = $oarBase."-O $logDir/bam2gvcfGatk.$sampsString.out -E $logDir/bam2gvcfGatk.$sampsString.err ";
	$oar .= "\" perl $bam2gvcf --indir $inDir --genome $genome --chroms $chroms --outdir $outDir --samples ";
	$oar .= join(',', @samples);
	$oar .= " --tmpdir /var/tmp/NTM_$sampsString --jobs $jobs --gatk \'$gatk\' --real\"";
	#warn "$oar\n";
	system($oar);
	@samples = ();
    }
}
