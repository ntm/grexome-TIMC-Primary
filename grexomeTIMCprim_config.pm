# NTM
# 18/09/2020

# Define the paths and filenames that should be install-specific and
# are needed by grexome-TIMC-primary.pl .
# Every hard-coded path/filename in the pipeline should be here.
# Also define some data-specific or behavioral config, if any.
#
# I am starting from grexomeTIMCsec_config.pm, several exported
# functions should be identical and in fact if you use both 
# primary and secondary pipelines you could have a single
# shared config.pm.


package grexomeTIMCprim_config;

use strict;
use warnings;
use Exporter;
our @ISA = ('Exporter');
# NOTE: this file can and should be copied somewhere and customized, 
# therefore we should never "use grexomeTIMCprim_config" but instead
# provide the customized *config.pm as an argument, see --config in
# grexome-TIMC-secondary.pl for an example.
our @EXPORT_OK = qw(refGenome refGenomeChromsBed fastTmpPath);


#################################################################
# files and paths needed by the pipeline 


# Return the reference human genome fasta file, with path
sub refGenome {
    # return the first file that exists, so this works on
    # both luxor and fauve
    foreach my $genome ("/home/nthierry/HumanGenome/hs38DH.fa",
			"/data/HumanGenome/hs38DH.fa",
			"/bettik/nthierry/HumanGenome/hs38DH.fa") {
	(-f $genome) && return($genome);
    }
    # if we get here no file was found...
    die "E: no refGenome found, you need to edit *config.pm";
}

# Full path to gzipped hs38_chrom-only BED file, with the full
# chromosomes 1-22,X,Y,M.
# This allows to ignore decoy/unplaced/alt regions. 
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
### tmpfs /mnt/RamDisk tmpfs size=192g 0 0
sub fastTmpPath {
    foreach my $ramdisk ("/mnt/RamDisk/", "/var/tmp") {
	(-d $ramdisk) && return($ramdisk);
    }
    # if we get here no dir was found...
    die "E: no fastTmpPath found, you need to edit *config.pm";
}


# module loaded ok
1;
