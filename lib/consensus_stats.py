#!/usr/bin/python
############################################################
# Generates base counts for consensus sequences
############################################################

import sys, os, re

############################################################

def readFasta(fa_file):
# A function to read a fasta file into a dictionary with each header as a key
    seqdict, titles = {}, [];
    for line in open(fa_file):
        if line == "\n":
            continue;
        line = line.strip();
        if line.startswith(">"):
            curkey = line;
            titles.append(curkey);
            seqdict[curkey] = "";
        else:
            seqdict[curkey] += line;
    return seqdict, titles;

############################################################
# Main

consensus_file, mode, output_file = sys.argv[1:];
# Parse positional input arguments:
# consensus_file: The combined consensus file for a given variant caller
# mode: The variant caller used: ivar or gatk
# output_file: The csv file to output the summary table with consensus stats

out_headers = "sample,A,T,C,G,N,length,perc.n,n.filter";
# The headers for the output file

with open(output_file, "w") as outfile:
    outfile.write(out_headers + "\n");
    # Open the output file and write the headers

    cur_seqs, cur_titles = readFasta(consensus_file);
    # Read the sequences

    for header in cur_seqs:
    # Go through every sequence

        if mode == "ivar":
            sample = header.split("_")[1];
        elif mode == "gatk":
            sample = header.split("_")[0];
        sample = sample.replace(">", "");
        # Get the sample ID based on the input mode

        a = cur_seqs[header].count("A") + cur_seqs[header].count("a");
        t = cur_seqs[header].count("T") + cur_seqs[header].count("t");
        c = cur_seqs[header].count("C") + cur_seqs[header].count("c");
        g = cur_seqs[header].count("G") + cur_seqs[header].count("g");
        n = cur_seqs[header].count("N") + cur_seqs[header].count("n");
        # Count each type of base in the current sequence

        seqlen = len(cur_seqs[header]);
        # Get the total sequence length

        # print(len(cur_seqs[header]))
        # print(sum([a,t,c,g,n]))

        # if "R" in cur_seqs[header]:
        #     print(cur_seqs[header].index("R"))

        # replacements = [(b, "") for b in "ATCGNatcgn"];
        # tmp = cur_seqs[header];
        # for pat, repl in replacements:
        #     tmp = re.sub(pat, repl, tmp);

        # print(header);
        # print(tmp);
        # print("----")
        # Some debug code...

        if mode == "gatk":
            assert seqlen == sum([a,t,c,g,n]), " * ERROR: Non-standard bases present: " + consensus_file + " " + header;
        # For sequences, make sure the sum of all base types is the sequence length
        # Don't do this for ivar because I've noticed that they can have ambiguous bases

        # Edge case: a sequencing failure that leads to a 0-length conesnsus
        # Might want to deal with this further upstream, but for now this fixes it
        if seqlen == 0:
            perc_n = "N/A";
        else:
            perc_n = n / seqlen;
        # Calculate the percentage of Ns

        if n > 5000:
            n_filter = "FILTER";
        else:
            n_filter = "PASS";
        # Determine whether this sequence passes the N filter

        outline = [ str(val) for val in [sample, a, t, c, g, n, seqlen, perc_n, n_filter] ];
        outfile.write(",".join(outline) + "\n");
        # Compile and write the output line

############################################################