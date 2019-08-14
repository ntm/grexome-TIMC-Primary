#!/usr/bin/perl


# 24/06/2019
# NTM

# Takes 2 args: inDir, outDir.
# grab the genome.S1.vcf.gz and associated .tbi files in each 
# grexomeXXXX subdir of $inDir, and move them to $outDir, 
# renaming to grexomeXXXX.g.vcf.gz .

# inDir: dir containing grexomeXXXX subdirs as produced by strelka
my $inDir;

# outDir: dir where GVCFs are moved
my $outDir;

(@ARGV == 2) || die "need 2 args: inDir outDir\n";

($inDir,$outDir) = @ARGV;

(-d $outDir) || (mkdir($outDir)) || 
    die "outDir $outDir doesn't exist and can't be mkdir'd\n";

opendir(IN, $inDir) || die "cannot opendir inDir $inDir\n";

while (my $grex = readdir(IN)) {
    ($grex =~ /^(grexome\d\d\d\d)$/) ||
	((warn "W: grex dir $grex doesn't look like a grexomedir, skipping\n") && next);
    my $inFile = "$inDir/$grex/results/variants/genome.S1.vcf.gz";
    foreach my $ext ("", ".tbi") {
	(-e "$inFile$ext") ||
	    ((warn "W: skipping $grex$ext because I can't find GVCF$ext file $inFile$ext\n") && next);
	my $outFile = "$outDir/$grex.g.vcf.gz$ext";
	(-e $outFile) &&
	    die "E: about to move GVCF$ext for $grex but it already exists in outDir $outDir, WTF!!\n";

	rename("$inFile$ext",$outFile) ||
	    die "E: rename failed for $inFile $outFile\n";
    }
}

closedir(IN);
