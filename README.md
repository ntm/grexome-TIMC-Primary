## Exome pipeline from TIMC in Grenoble - primary analyses: from FASTQs to a single multi-sample GVCF per variant-caller.

The pipeline currently performs variant-calling independently with Strelka2, GATK4, elPrep5, and deepVariant.


*****************
### INSTALLATION:
```
git clone https://github.com/ntm/grexome-TIMC-Primary.git
```

*****************
### DEPENDENCIES:
We try to keep external dependencies to a minimum.

#### PERL modules
The only required PERL modules are listed below. Most are standard core modules and should already be available on your system, a few will probably need to be installed from your distribution's standard repositories (e.g. "sudo yum install perl-Parallel-ForkManager perl-Spreadsheet-XLSX" on RHEL7).
  - Standard core modules:
    - Exporter
    - File::Basename
    - File::Copy
    - File::Glob
    - File::Path
    - File::Spec
    - File::Temp
    - FindBin
    - Getopt::Long
    - POSIX
    - Cwd
  - Other modules:
    - Parallel::ForkManager
    - Spreadsheet::XLSX

#### External tools
We use the following external tools. They must be installed and should be in your PATH. If they are not in PATH, you can customize them by providing additional options to each script, eg `1_fastq2bam.pl --binpath`. This will require editing the master script `grexome-TIMC-primary.pl`, you can open an issue if you need help.
- grep, zgrep, gunzip
- bgzip, tabix [from HTSlib](http://www.htslib.org/download/)
- [fastp](https://github.com/OpenGene/fastp)
- [bwa](https://github.com/lh3/bwa)  OR [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
- [samblaster](https://github.com/GregoryFaust/samblaster)
- [samtools](http://www.htslib.org/download/)
- bwa-postalt.js and k8 from Heng Li's [BWA-kit package](https://sourceforge.net/projects/bio-bwa/files/bwakit/)

In addition you may need the following, depending on the variant-callers you specify with --callers:
- configureStrelkaGermlineWorkflow.py from [STRELKA2](https://github.com/Illumina/strelka), for 2_bam2gvcf_strelka.pl
- gatk from [GATK4](https://github.com/broadinstitute/gatk/), for 2_bam2gvcf_gatk.pl
- elprep from [ELPREP5](https://github.com/exascience/elprep), for 2_bam2gvcf_elprep.pl
- singularity and a [deepvariant](https://github.com/google/deepvariant) singularity image (SIF format), for 2_bam2gvcf_deepVariant.pl

Refer to each tool's github / webpage instructions for installation.


*****************
### REQUIRED DATA:
The reference genome in fasta format is required, and it must be indexed/preprocessed as needed by bwa, bwa-mem2, strelka, gatk, elprep and deepVariant. Suggested protocol:
- download the latest version of [BWA-kit](https://sourceforge.net/projects/bio-bwa/files/bwakit/)
- edit run-gen-ref to [update the url38= line](https://github.com/lh3/bwa/issues/189)
- produce hs38DH.fa and index/preprocess it as follows:
```
run-gen-ref hs38DH
bwa index hs38DH.fa ## for BWA
bwa-mem2 index hs38DH.fa ## for BWA-MEM2 if you use it, careful this indexing requires tons of RAM
samtools faidx hs38DH.fa ## for strelka (and deepVariant?)
gatk CreateSequenceDictionary -R hs38DH.fa ## for GATK
elprep fasta-to-elfasta hs38DH.fa hs38DH.elfasta --log-path . ## for elprep
```

In addition, a small BED file with chromosomes 1-22,X,Y,M is required for the variant callers.

To produce this file, you can copy-paste [this content](https://github.com/Illumina/strelka/blob/v2.9.x/docs/userGuide/README.md#improving-runtime-for-references-with-many-short-contigs-such-as-grch38) into a file (eg hs38_chroms.bed), then bgzip and index it:
```
bgzip hs38_chroms.bed
tabix -p bed hs38_chroms.bed.gz
```
Alternately a copy of the resulting bgzipped file is [provided](Metadata/hs38_chroms.bed.gz) (still needs to be tabix-indexed).


*****************
### CONFIGURATION:
Before using the pipeline you must customize (a copy of) the grexomeTIMCprim_config.pm file, which defines every install-specific variable (eg path+name of the reference human genome fasta file, possibly produced following the instructions in "REQUIRED DATA" above).

Every subroutine in grexomeTIMCprim_config.pm is self-documented and will need to be customized.

To do this you should copy the file somewhere and edit the copy, then use `--config`. Otherwise you could just edit the distributed copy in-place, although this not as flexible and your customizations will cause conflicts when you git pull.


******************
### METADATA FILES:
The pipeline uses a single metadata XLSX file, which describes the samples. This metadata file (and the code that parses it, in grexome_metaParse.pm) is shared with the [grexome-TIMC-Secondary pipeline](https://github.com/ntm/grexome-TIMC-Secondary). Therefore some columns are required in the XLSX despite not being used in this primary pipeline. The provided [example file](Metadata/samples.xlsx) can serve as a starting point. All columns present in this file (and listed below) are required, don't change their names! You can add new columns and/or change the order of columns to your taste, just don't touch the pre-existing column names. 

Required columns:
- sampleID: unique identifier for each sample, these are typically created with a uniform naming scheme when new samples are integrated into the pipeline and are used internally throughout the pipeline. Rows with sampleID=="0" are ignored, this allows to retain metadata info about obsolete samples.
- specimenID: external identifier for each sample, typically related to the original FASTQ filenames.
- patientID: can be empty, ignored in grexome-TIMC-Primary.
- pathologyID: required (alphanumeric string) but ignored in grexome-TIMC-Primary.
- Causal gene: can be empty, ignored in grexome-TIMC-Primary.

Optional column:
- Sex: if this column exists each sample must be 'F' or 'M', used by 0_qc_checkSexChroms.pl .



*****************
### EXAMPLE USAGE:
```
DATE="220525"
perl ~/grexome-TIMC-Primary/grexome-TIMC-primary.pl --samples=samples.xlsx --workdir=PrimaryAnalyses_${DATE} --callers=strelka,gatk,elprep --config=myPrimConfig.pm &> grexomeTIMCprim_${DATE}.log &
```

This will "process" every sampleID in samples.xlsx, ie:
- produce analysis-ready BAMs (trim, align, mark dupes, sort);
- produce individual GVCFs with STRELKA2, GATK4 and ELPREP5;
- fix variant-caller quirks and bugs, filter low-quality variant calls;
- produce a single multi-sample merged GVCF per variant-caller;
- produce QC files analyzing HOMO/HET calls on the sex chromosomes (optional, only if samples.xlsx has a "Sex" column).

For each sample, any step where the result file already exists is skipped.

FASTQs, BAMs and GVCFs are searched for / produced in a hierarchy of subdirs defined at the top of grexome-TIMC-primary.pl, all within \$dataDir (defined in your customized copy of grexomeTIMCprim_config.pm, see CONFIGURATION).

Logs and copies of the metadata are produced in the provided workdir (PrimaryAnalyses_${DATE}/).

More extensive documentation is obtained with:
```
perl grexome-TIMC-primary.pl --help
```


**********************
### OTHER REPO CONTENT:

##### PREP STEP #####
0_grexomizeFastqs.pl : house-keeping script used to organize our FASTQ files before running the grexome-TIMC-primary pipeline. It parses the samples.xlsx metadata file and takes original fastqs to produce a consolidated uniform dataset of "grexomized" fastqs. It's probably not re-usable as-is since it depends on your file-naming conventions, but should provide a good starting point.


##### RUNNING SOME STEPS ON CLUSTERS #####
- 1_fastq2bam_makeOarJobs.pl and README_fastq2bam_makeOarJobs :
- Bam2gvcf_Strelka_PackagedWithBinaries/ :
- Bam2gvcf_GATK_PackagedWithBinaries/ :

scripts + READMEs for running fastq2bam, bam2gvcf_strelka or bam2gvcf_gatk on other systems / compute clusters, eg f-dahu. Probably not re-usable as-is, but may be useful.


##### EXPERIMENTAL STEPS, NOT USED #####
- 4_CombineGVCFs_GATK/ : attempt to use GATK CombineGVCFs, crazy slow.
- 4_GVCFs2GenomicsDB/ : import single-sample GATK GVCFs into a GenomicsDB.
- 5_GenomicsDB2VCF/ : perform joint genotyping with GATK GenotypeGVCFs from a GenomicsDB. Unusably slow.

