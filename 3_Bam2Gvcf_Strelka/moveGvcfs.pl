#!/usr/bin/perl


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
	    die "E $0: rename failed for $inFile $outFile\n";
    }
}

closedir(IN);
