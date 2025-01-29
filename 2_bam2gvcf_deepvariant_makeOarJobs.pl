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


# 14/06/2022
# NTM
#
# submit bam2gvcf_deepvariant.pl jobs on dahu with oar.
# Not sure if this script can be useful outside our organization:
# it would need adapting for a different batch scheduler and/or
# cluster, paths and filename patterns are hard-coded (eg samples
# are called grexome\d\d\d\d), etc... 
#
# deepvariant is executed via singularity.
#
# Takes 2 args: $first and $last, the first and last grexomes (samples)
# for this run of bam2gvcf_gatk_makeOarJobs.pl
# NOTE: not more than 50*$grexomesPerJob at a time, or you get error:
# Admission Rule ERROR : [ADMISSION RULE] Error: you cannot have more than 50 jobs waiting in the queue at the same time.


use strict;
use warnings;

# max number of grexomes to process in a single OAR job of 24h, 
# first tests on luxor with 30 jobs took ~8h10 to process 6 samples, dahu nodes
# should be a bit faster (newer CPUs) but some samples may take longer,
# to play safe let's limit at 14
my $grexomesPerJob = 14;

(@ARGV == 2) || die "E: need two ints as arguments, first and last\n";
my ($first,$last) = @ARGV;

($last >= $first) || 
    die "E: must have first <= last, here we have: first $first last $last\n";

(($last - $first) > 50*$grexomesPerJob) &&
    die "E: cannot queue more than 50 jobs at a time on OAR/dahu\n";


#################################
# hard-coded stuff

# number of parallel threads for one sample on one node
my $jobs = 32;

# hard-coded dirs:
# log to $logDir
my $logDir = "/bettik/thierryn/ProcessBams_DV_2206_logs/";
# bams are in $inDir
my $inDir = "/bettik/thierryn/BAMs_All_Selected/";
# produce gvcfs in $outDir
my $outDir = "/bettik/thierryn/ProcessBams_DV_2206/";
# scratch on f-dahu nodes: /var/tmp
my $tmpDir = "/var/tmp/DVtmp/";
# path to deepvariant SIF
my $deepvariant = "/home/thierryn/Software/DeepVariant/deepvariant_1.4.0.sif";

# path to bam2gvcf_deepvariant.pl
my $bam2gvcf = "/home/thierryn/Software/grexome-TIMC-Primary/2_bam2gvcf_deepvariant.pl";
# ref genome
my $genome = "/bettik/thierryn/HumanGenome/hs38DH.fa";
# chroms or target regions
my $chroms = "/bettik/thierryn/HumanGenome/hs38_chroms.bed.gz";

# oarsub command with params: run on my project ngs-timc, and...
my $oarBase = "oarsub --project ngs-timc";
## INITIAL BIG JOBS TO PROCESS GREXOMES 50-630: ask for 1 full node, 24h walltime, per job
$oarBase .= " -l /nodes=1,core=32,walltime=24 ";

my $grex = $first;
while ($grex <= $last) {
    my $slotsThisJob = $grexomesPerJob;
    my $grexStartJob = $grex;
    # build comma-separated list of samples
    my $samples = "";
    while ($slotsThisJob > 0) {
	($grex <= $last) || last;
	my $s = $grex;
	$grex++;
	# left-pad with zeroes to 4 digits
	($s < 10) && ($s = "0$s");
	($s < 100) && ($s = "0$s");
	($s < 1000) && ($s = "0$s");
	# skip if BAM doesn't exist
	(-e "$inDir/grexome$s.bam") ||
	    ((warn "W: $inDir/grexome$s.bam doesn't exist, skipping $s\n") && next);
	# ok add this sample and decrement remaining slots for this job
	$samples .= "grexome$s,";
	$slotsThisJob--;
    }
    # chop trailing comma
    (chop($samples) eq ',') ||
	((warn "W: no remaining samples starting at $grexStartJob, all done\n") && last);
    my $oar = $oarBase." -O $logDir/bam2gvcf.$grexStartJob.out -E $logDir/bam2gvcf.$grexStartJob.err ";
    $oar .= "\" perl $bam2gvcf --indir $inDir --samples $samples --genome $genome --chroms $chroms ";
    $oar .= "--outdir $outDir --tmpdir $tmpDir --deepvariant $deepvariant --jobs $jobs --real\"";


    #warn "$oar\n";
    system($oar);
}
