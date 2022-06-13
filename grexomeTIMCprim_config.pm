# NTM
# 18/09/2020

# Define the paths and filenames that should be install-specific and
# are needed by grexome-TIMC-primary.pl .
# Every hard-coded path/filename in the pipeline should be here.
# Also define some data-specific or behavioral config, if any.


package grexomeTIMCprim_config;

use strict;
use warnings;
use Exporter;
our @ISA = ('Exporter');
# NOTE: this file can be copied somewhere and customized, therefore
# we never "use grexomeTIMCprim_config" but instead provide the
# customized *config.pm as an argument, see --config in grexome-TIMC-primary.pl
# for an example.
our @EXPORT_OK = qw(dataDir fastqDir mirror refGenome refGenomeElPrep refGenomeChromsBed 
		    fastTmpPath binPath bwakitPath deepVariantSIF);


#################################################################
# files and paths needed by the pipeline 

# dir holding the hierarachy of subdirs and files that will be
# populated with all results (BAMs, GVCFs).
# The hierarchy (hard-coded in grexome-TIMC-primary.pl) doesn't need
# to change, but &dataDir() certainly does.
sub dataDir {
    my $dataDir = "/data/nthierry/PierreRay/";
    return($dataDir);
}

# dir containing the "grexomized" FASTQs: for each sample we expect
# a single pair of FASTQ files in $fastqDir, and these files must be
# named ${sample}_1.fq.gz and ${sample}_2.fq.gz .
# If you have several pairs of files for some samples and/or heterogeneous filename
# conventions, you need to preprocess your FASTQs (merge and/or rename and/or symlink).
# We use 0_grexomizeFastqs.pl for this, it could be helpful but you'll likely have to
# adapt it so it meets your local conventions.
sub fastqDir {
    #default to a subdir of &dataDir()
    my $fastqDir = &dataDir()."/FASTQs_All_Grexomized/";
    return($fastqDir);
}


# rsync path where you maintain a mirror of the BAMs and GVCFs.
# This string is never executed, it only appears at the end of the logs to
# allow easy copy-pasting.
# Return "" to NOT print these final log lines (eg if you mirror via cron)
sub mirror {
    my $mirror = "cargo.u-ga.fr:/bettik/thierryn/";
    return($mirror);
}


# Return the reference human genome fasta file, with path.
# It needs to be indexed/preprocessed in various ways for all
# the tools we use, see "REQUIRED DATA" section in the README.
sub refGenome {
    # return the first file that exists, so this works on
    # all our servers
    foreach my $genome ("/home/nthierry/HumanGenome/hs38DH.fa",
			"/data/HumanGenome/hs38DH.fa",
			"/bettik/nthierry/HumanGenome/hs38DH.fa") {
	(-f $genome) && return($genome);
    }
    # if we get here no file was found...
    die "E: no refGenome found, you need to edit *config.pm";
}

# Return the reference human genome in the elPrep5 elfasta format, with path
# This can be produced with eg:
# elprep fasta-to-elfasta hs38DH.fa hs38DH.elfasta --log-path .
sub refGenomeElPrep {
    # return the first file that exists, so this works on
    # all our servers
    foreach my $genome ("/home/nthierry/HumanGenome/hs38DH.elfasta",
			"/data/HumanGenome/hs38DH.elfasta",
			"/bettik/nthierry/HumanGenome/hs38DH.elfasta") {
	(-f $genome) && return($genome);
    }
    # if we get here no file was found...
    die "E: no refGenome.elfasta found, you need to edit *config.pm";
}

# Full path to bgzipped and tabix-indexed hs38_chrom-only BED file,
# with the full chromosomes 1-22,X,Y,M.
# This allows to ignore decoy/unplaced/alt regions.
# An example of the expected file is provided in Metadata/ .
sub refGenomeChromsBed {
    # default to a file alongside the ref genome
    my $chromsBed = &refGenome();
    ($chromsBed =~ s/hs38DH.fa$/hs38_chroms.bed.gz/) ||
	die "E: cannot substitute hs38 fasta for bed filename, maybe refGenome isn't named as expected?";
    (-f $chromsBed) ||
	die "E: chromsBed file $chromsBed doesn't exist, rsync it from somewhere or fix the code";
    return($chromsBed);
}


# Return a tmp dir with fast RW access, ideally a ramdisk, otherwise
# hopefully at least an SSD.
# Example: on linux you can make a 96GB ramdisk (assuming you have enough RAM)
# accessible in /mnt/RamDisk/ with:
### sudo mkdir /mnt/RamDisk
### sudo mount -t tmpfs -o size=96g myramdisk /mnt/RamDisk/
# You can make the mount automatic on boot by adding to /etc/fstab:
### tmpfs /mnt/RamDisk tmpfs size=96g 0 0
sub fastTmpPath {
    foreach my $ramdisk ("/mnt/RamDisk/", "/var/tmp") {
	(-d $ramdisk) && return($ramdisk);
    }
    # if we get here no dir was found...
    die "E: no fastTmpPath found, you need to edit *config.pm";
}


# Path to programs used (fastp, samblaster, samtools, bwa-mem2 / bwa, ),
# Returning the empty string searches in $PATH.
# If your binaries are in various dirs all over the place, you should
# create a subdir somewhere and symlink all the required binaries there,
# then use that subdir as binPath
sub binPath {
    # default to searching in PATH
    my $binPath = "";
    return($binPath);
}

# Dir containing Heng Li's bwa-kit package - we use k8 and bwa-postalt.js to correctly
# deal with ALT contigs.
# See https://github.com/lh3/bwa/blob/master/README-alt.md
sub bwakitPath {
    # path to bwa-kit subdir (with k8 and bwa-postalt.js), available on SF:
    # https://sourceforge.net/projects/bio-bwa/files/bwakit/
    my $bwakit = "/home/nthierry/Software/BWA-kit/bwa.kit/";
    (-d $bwakit) ||
	die "E: specified bwakitPath $bwakit not found, you need to edit *config.pm";
    return($bwakit);
}


# Return the path/name of the deepVariant singularity image you want to use (assuming
# you want to include deepVariant in --callers).
# For example, we currently produce our image on centos7 with:
### BIN_VERSION="1.4.0"
### singularity pull docker://google/deepvariant:"${BIN_VERSION}"
sub deepVariantSIF {
    my $dvSif = "/home/nthierry/Software/DeepVariant/deepvariant_1.4.0.sif";
    (-f $dvSif) ||
	die "E: deepVariant singularity image $dvSif not found, you need to edit *config.pm";
    return($dvSif);
}


# module loaded ok
1;
