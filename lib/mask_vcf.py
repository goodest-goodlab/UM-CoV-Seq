#!/usr/bin/python
############################################################
# Masks problematic SARS-CoV-2 sites in vcf files by
# marking them with the FILTER tag.
############################################################

import sys, re, os

############################################################

def readVCF(vcffile):
# Reads a vcf file .
    iupac = {"A":["A"], "T":["T"], "C":["C"], "G":["G"], "R":["A","G"], "Y":["C","T"], "S":["G","C"], "W":["A","T"], "K":["G","T"], "M":["A","C"], 
                    "B":["C","G","T"], "D":["A","G","T"], "H":["A","C","T"], "V":["A","C","G"], "N":["A","T","C","G"], ".":["-"] };

    vcflines = [ line.strip().split("\t") for line in open(vcffile) if not line.startswith("#") ];
    for i in range(len(vcflines)):
        nts = vcflines[i][3].split(",");
        new_nts = [];
        for nt in nts:
            new_nts += iupac[nt];
        vcflines[i][3] = list(set(new_nts));

    return vcflines;

############################################################

input_file, output_file = sys.argv[1:];

problematic_vcf = "ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf";
problematic_sites = readVCF(problematic_vcf);

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

    for site in problematic_sites:
        if vcflines[i][1] == site[1] and vcflines[i][4] in site[4] and vcflines[i][6] == "PASS":
            vcflines[i][6] = "FILTER";
            num_filtered += 1;
    # Check each SNP in the provided -vcf file. If it matches the current SNP, add the filter string to the
    # FILTER column.

with open(output_file_unzipped, "w") as new_vcf:
    for line in vcflines:
        new_vcf.write("\t".join(line) + "\n");
# Re-write the iteration VCF file.

os.system("bgzip " + input_file_unzipped);
os.system("bgzip " + output_file_unzipped);
# Re-compress the iteration VCF files.

############################################################