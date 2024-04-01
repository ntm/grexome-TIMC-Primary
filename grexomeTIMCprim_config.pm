
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
		    fastTmpPath binPath bwakitPath strelkaBin deepVariantSIF gatkBin elprepBin);


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
    #my $mirror = "cargo.u-ga.fr:/bettik/thierryn/";
    my $mirror = "";
    return($mirror);
}


# Return the reference genome fasta file, with path.
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

# Return the reference genome in the elPrep5 elfasta format, with path
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

# Full path to bgzipped and tabix-indexed BED file specifying regions where variants should
# be called, other regions are ignored.
# For example this can be a hs38_chrom-only BED file with the full chromosomes 1-22,X,Y,M:
# this allows to ignore decoy/unplaced/alt regions, one such file is provided in Metadata/ .
# Return "" to call variants everywhere.
sub refGenomeChromsBed {
    # default: call variants on the whole genome
    my $chromsBed = "";
    # below is what we use for GRCh38, expecting a BED file with a specific name alongside the ref genome...
    # It most likely needs to be adapted to your file naming conventions, or just commented out
    my $chroms = &refGenome();
    if ($chroms !~ s/hs38DH.fa$/hs38_chroms.bed.gz/) {
	warn "W: in refGenomeChromsBed, cannot build BED filename, you should edit grexomeTIMCprim_config.pm";
    }
    elsif (! -f $chroms) {
	warn "W: in refGenomeChromsBed, $chroms doesn't exist, you should edit grexomeTIMCprim_config.pm";
    }
    else {
	$chromsBed = $chroms;
    }
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

# Dir containing Heng Li's bwa-kit package - when aligning on GRCh38, we use k8 and
# bwa-postalt.js to correctly deal with ALT contigs.
# See https://github.com/lh3/bwa/blob/master/README-alt.md
# Return '' to skip the bwakit-postAlt step (eg when working with non-human data)
sub bwakitPath {
    # path to bwa-kit subdir (with k8 and bwa-postalt.js), available on SF:
    # https://sourceforge.net/projects/bio-bwa/files/bwakit/
    my $bwakit = "/home/nthierry/Software/BWA-kit/bwa.kit/";
    (-d $bwakit) ||
	die "E: specified bwakitPath $bwakit not found, you need to edit *config.pm";
    return($bwakit);
}


# path+name of configureStrelkaGermlineWorkflow.py from strelka distrib,
# or use a two-liner wrapper that can be in your PATH or in &binPath(),
# such as the following (which we have in /usr/local/bin/strelkaGermline.sh):
### #!/bin/sh
### /home/nthierry/Software/Strelka/strelka-latest/bin/configureStrelkaGermlineWorkflow.py "$@"
sub strelkaBin {
    my $strelka = "strelkaGermline.sh";
    # prepend &binPath() if it's non-empty and $strelka exists there
    (&binPath() ne "") && (-f &binPath()."/$strelka") && ($strelka = &binPath()."/$strelka");
    return($strelka);
}

# path+name of the deepVariant singularity image you want to use (assuming
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

# path+name of GATK launcher script distributed with GATK4, default is "gatk"
sub gatkBin {
    my $gatk = "gatk";
    # prepend &binPath() if it's non-empty and $gatk exists there
    (&binPath() ne "") && (-f &binPath()."/$gatk") && ($gatk = &binPath()."/$gatk");
    return($gatk);
}

# path+name of elPrep5 binary, default is "elprep"
sub elprepBin {
    my $elprep = "elprep";
    # prepend &binPath() if it's non-empty and $elprep exists there
    (&binPath() ne "") && (-f &binPath()."/$elprep") && ($elprep = &binPath()."/$elprep");
    return($elprep);
}


# module loaded ok
1;
