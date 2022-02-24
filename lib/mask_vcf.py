#!/usr/bin/python
############################################################
# Masks problematic SARS-CoV-2 sites in vcf files by
# marking them with the FILTER tag.
############################################################

import sys, re, os

############################################################

def readVCF(vcffile):
# Reads a vcf file.
    iupac = {"A":["A"], "T":["T"], "C":["C"], "G":["G"], "R":["A","G"], "Y":["C","T"], "S":["G","C"], "W":["A","T"], "K":["G","T"], "M":["A","C"], 
                    "B":["C","G","T"], "D":["A","G","T"], "H":["A","C","T"], "V":["A","C","G"], "N":["A","T","C","G"], ".":["-"] };
    # The problematic sites are defined with IUPAC codes, so we include them here

    vcflines = [ line.strip().split("\t") for line in open(vcffile) if not line.startswith("#") ];
    for i in range(len(vcflines)):
        nts = vcflines[i][3].split(",");
        new_nts = [];
        for nt in nts:
            new_nts += iupac[nt];
        vcflines[i][3] = list(set(new_nts));
    # This reads the lines of the vcf file, replacing any ambiguities in the alternate allele with
    # a list of nucleotides.

    return vcflines;

############################################################

input_file, output_file = sys.argv[1:];
# Parse positional input arguments:
# input_file: An input VCF file
# output_file: An output, masked VCF file

input_file_unzipped = input_file.replace(".gz", "");
output_file_unzipped = output_file.replace(".gz", "");
os.system("gunzip " + input_file);
vcflines = [ line.strip().split("\t") for line in open(input_file_unzipped) ];
# Unzip and read the iteration VCF file.

num_filtered = 0;
for i in range(len(vcflines)):
    if vcflines[i][0].startswith("#"):
        continue;
    # Check each SNP in the VCF file; skip the header lines.

    genotype = vcflines[i][9].split(":")[0];
    if genotype == ".":
        vcflines[i][4] = "N";
        vcflines[i][6] = "PASS";
    # To mask sites without enough info for a call, switch their alt allele to N and the filter to PASS to ensure
    # the N is inserted into the consensus.

with open(output_file_unzipped, "w") as new_vcf:
    for line in vcflines:
        new_vcf.write("\t".join(line) + "\n");
# Re-write the iteration VCF file.

os.system("bgzip " + input_file_unzipped);
os.system("bgzip " + output_file_unzipped);
# Re-compress the iteration VCF files.

############################################################
