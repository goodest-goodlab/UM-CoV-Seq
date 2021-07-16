#!/usr/bin/python
############################################################
# Masks problematic SARS-CoV-2 sites in pileup files by
# changing all mapped bases at those sites to the reference
# to prevent variants from being called.
############################################################

import sys, re

############################################################

# def readVCF(vcffile):
# # Reads a vcf file .
#     iupac = {"A":["A"], "T":["T"], "C":["C"], "G":["G"], "R":["A","G"], "Y":["C","T"], "S":["G","C"], "W":["A","T"], "K":["G","T"], "M":["A","C"], 
#                     "B":["C","G","T"], "D":["A","G","T"], "H":["A","C","T"], "V":["A","C","G"], "N":["A","T","C","G"], ".":["-"] };

#     sites = [ line.strip().split("\t")[1] for line in open(vcffile) if not line.startswith("#") ];

#     vcflines = [ line.strip().split("\t") for line in open(vcffile) if not line.startswith("#") ];
#     for i in range(len(vcflines)):
#         nts = vcflines[i][3].split(",");
#         new_nts = [];
#         for nt in nts:
#             new_nts += iupac[nt];
#         vcflines[i][3] = list(set(new_nts));

#     return vcflines;

############################################################

input_file, output_file = sys.argv[1:];

problematic_vcf = "ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf";
problematic_sites = [ line.strip().split("\t")[1] for line in open(problematic_vcf) if not line.startswith("#") ];

replacements = [(b, ".") for b in "ATCGNatcgn"];

with open(output_file, "w") as outfile:
    for line in open(input_file):
        line_list = line.strip().split("\t");
        if line[1] not in problematic_sites:
            outfile.write(line);
        else:
            #ref = line_list[2];
            base_str = line_list[4];
            for pat, repl in replacements:
                base_str = re.sub(pat, repl, base_str);

            line_list[4] = base_str;

            outfile.write("\t".join(line_list) + "\n");

############################################################