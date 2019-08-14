#!/usr/bin/perl


# 17/05/2019
# NTM

# Takes 2 args: a subdir containing input GVCFs ending in g.vcf.gz, 
# and an $outFile which will be produced: a single vcf.gz file.
# Any GATK errors/warnings/logging is printed to stderr.

use strict;
use warnings;
use Parallel::ForkManager;

my $numJobs = 12 ;
my $pm = new Parallel::ForkManager($numJobs);







# Take one BAM file (must exist) and the name of a GVCF file (will be
# squashed if it exists).
# Produce the GVCF with GATK.
sub bam2gvcf_GATK {
    # absolute path to GATK .jar (need abs path so can be called anywhere)
    my $gatk = "/home/nthierry/Software/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar";

# absolute path to reference genome fasta
my $refGenome = "/home/nthierry/Grexome/ressources/GRCh38/grch38_essential.fa";

# absolute path to tmpdir with a lot of free space
my $tmpDir = "/home/nthierry/VariantCalling/tmp/" ;

# max memory given to java (end with G for GBs)
my $maxMem = "96G";

# number of GATK threads
my $numThreads = 12;


(@ARGV == 2) || 
    die "needs 2 args: a subdir containing input GVCFs and an outFile to produce\n";
my ($inDir, $outFile) = @ARGV;

(-d $inDir) || 
    die "inDir $inDir is not a directory\n";
(-f $outFile) && 
    die "outFile $outFile already exists, remove it or choose another outFile name\n";
($outFile =~ /\.vcf\.gz$/) ||
    die "outFile $outFile should end in .vcf.gz or GATK probably won't work / gzip it.\n";

opendir(INDIR, $inDir) ||
    die "cannot opendir inDir $inDir\n";


example usage:
 java -jar GenomeAnalysisTK.jar \
     -R reference.fasta \
     -T HaplotypeCaller \
     -I sample1.bam \
     --emitRefConfidence GVCF \
     -o output.raw.snps.indels.AS.g.vcf

 -nct broken for many users, avoid it

example from history on yeast:
 java -Xmx16G -jar /home/nthierry/Software/GATK/GenomeAnalysisTK-3.8-1-0-gf15c1c3ef/GenomeAnalysisTK.jar -T HaplotypeCaller -R ../RefGenome/GCF_000146045.2_R64_genomic.fna -I ../AlignFastqsAndSortBAMs/A.S288C.sorted.bam -ERC BP_RESOLUTION -o A.S288C.g.vcf -dt NONE -nct 8 -S STRICT -A StrandAlleleCountsBySample


more info:
https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_tools_walkers_haplotypecaller_HaplotypeCaller.php
https://software.broadinstitute.org/gatk/documentation/tooldocs/3.8-0/org_broadinstitute_gatk_engine_CommandLineGATK.php




# -Xms$maxMem -> starting mem to use, doesn't seem useful
my $command = "java -Xmx$maxMem -Djava.io.tmpdir=$tmpDir";
$command .= " -jar $gatk -nt $numThreads -S STRICT -T GenotypeGVCFs -R $refGenome";

foreach my $inFile (sort readdir(INDIR)) {
    ($inFile =~ (/\.g\.vcf\.gz$/)) || next;
    $command .= " -V $inDir/$inFile" ;
}

$command .= " -o $outFile";


system($command);
