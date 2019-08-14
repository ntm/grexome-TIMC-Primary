#!/usr/bin/perl


# 05/04/2018
# NTM

# Takes 2 args: a subdir containing input GVCFs .gz, and an $outFile
# which will be produced (merging all input GVCFs into a single GVCF).
# Any GATK errors/warnings/logging is printed to stderr.

# absolute path to GATK .jar (need abs path so can be called anywhere)
my $gatk = "/home/nthierry/Software/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";

# absolute path to reference genome fasta
my $refGenome = "/home/nthierry/Grexome/ressources/GRCh38/grch38_essential.fa";

(@ARGV == 2) || 
    die "needs 2 args: a subdir containing input GVCFs and an outFile to produce\n";
my ($inDir, $outFile) = @ARGV;

(-d $inDir) || 
    die "inDir $inDir is not a directory\n";
(-f $outFile) && 
    die "outFile $outFile already exists, remove it or choose another outFile name\n";

opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";



my $command = "java -Xmx96G -jar $gatk -T CombineGVCFs -R $refGenome";

while (my $inFile = readdir(INDIR)) {
    ($inFile =~ (/\.g\.vcf\.gz$/)) || next;
    $command .= " -V $inDir/$inFile" ;
}

$command .= " -o $outFile";


system($command);
