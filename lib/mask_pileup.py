#!/usr/bin/python
############################################################
# Masks problematic SARS-CoV-2 sites in pileup files by
# changing all mapped bases at those sites to the reference
# to prevent variants from being called.
############################################################

import sys, re

############################################################
# Maind

input_file, output_file = sys.argv[1:];
# Parse the positional input arguments:
# input_file: A pileup file
# output_file: A masked pileup file

problematic_vcf = "ProblematicSites_SARS-CoV2/problematic_sites_sarsCov2.vcf";
problematic_sites = [ line.strip().split("\t")[1] for line in open(problematic_vcf) if not line.startswith("#") ];
# Read the problematic sites from the specified VCF

replacements = [(b, ".") for b in "ATCGNatcgn"];
# Compile tuples of replacements
# For a pileup, we replace any non reference base at the problematic sites with the reference base, represented as "."

with open(output_file, "w") as outfile:
    for line in open(input_file):
        line_list = line.strip().split("\t");
        # Parse every line in the input pileup file

        if line_list[1] not in problematic_sites:
            outfile.write(line);
        # If the line is not one of the problematic sites, write it as-is to the output file
        else:
            base_str = line_list[4];
            # Get the base string

            for pat, repl in replacements:
                base_str = re.sub(pat, repl, base_str);
            # Replace all non reference bases in the base string by looping through all possibilities

            line_list[4] = base_str;
            # Replace the old base string with the new one

            outfile.write("\t".join(line_list) + "\n");
            # Write the modified line to the output pileup file
        # Replace alternate alleles at problematic sites

############################################################
