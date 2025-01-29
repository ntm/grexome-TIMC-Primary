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


# 24/06/2019
# NTM

# Takes 2 args: inDir, outDir.
# inDir contains one or more subdirs whose names are sampleIDs, each
# containing the results (files and subdirs) of a STRELKA run.
# 
# Grab the genome.S1.vcf.gz and associated .tbi files in each 
# $sample subdir of $inDir, and move them to $outDir, 
# renaming to $sample.g.vcf.gz .

use strict;
use warnings;
use File::Basename qw(basename);

# we use $0 in every stderr message but we really only want
# the program name, not the path
$0 = basename($0);


# inDir: dir containing $sampleID subdirs as produced by strelka
my $inDir;

# outDir: dir where GVCFs are moved
my $outDir;

(@ARGV == 2) || die "E $0: need 2 args: inDir outDir\n";

($inDir,$outDir) = @ARGV;

(-d $outDir) || (mkdir($outDir)) || 
    die "E $0: outDir $outDir doesn't exist and can't be mkdir'd\n";

opendir(IN, $inDir) || die "E $0: cannot opendir inDir $inDir\n";

while (my $sample = readdir(IN)) {
    # silently ignore anything starting with '.'
    ($sample =~ /^\./) && next;
    my $inFile = "$inDir/$sample/results/variants/genome.S1.vcf.gz";
    foreach my $ext ("", ".tbi") {
	(-e "$inFile$ext") ||
	    ((warn "W $0: skipping $sample$ext because I can't find GVCF$ext file $inFile$ext\n") && next);
	my $outFile = "$outDir/$sample.g.vcf.gz$ext";
	(-e $outFile) &&
	    die "E $0: about to move GVCF$ext for $sample but it already exists in outDir $outDir, WTF!!\n";

	rename("$inFile$ext",$outFile) ||
	    die "E $0: rename failed for $inFile$ext $outFile\n";
    }
}

closedir(IN);
