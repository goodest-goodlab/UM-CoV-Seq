#!/usr/bin/python
############################################################
# Generates base counts for consensus sequences
############################################################

import sys, os, re

############################################################

def readFasta(fa_file):
    seqdict = {};
    for line in open(fa_file):
        if line == "\n":
            continue;
        line = line.strip();
        if line.startswith(">"):
            curkey = line;
            seqdict[curkey] = "";
        else:
            seqdict[curkey] += line;
    return seqdict;

############################################################
# Main

consensus_file, mode, output_file = sys.argv[1:];
out_headers = "sample,seq.header,A,T,C,G,N,length,perc.n,n.filter";
with open(output_file, "w") as outfile:
    outfile.write(out_headers + "\n");

    cur_seqs = readFasta(consensus_file);
    for header in cur_seqs:

        if mode == "ivar":
            sample = header.split("_")[1];
        elif mode == "gatk":
            sample = header.split("_")[0];

        a = cur_seqs[header].count("A") + cur_seqs[header].count("a");
        t = cur_seqs[header].count("T") + cur_seqs[header].count("t");
        c = cur_seqs[header].count("C") + cur_seqs[header].count("c");
        g = cur_seqs[header].count("G") + cur_seqs[header].count("g");
        n = cur_seqs[header].count("N") + cur_seqs[header].count("n");

        seqlen = len(cur_seqs[header]);

        # print(len(cur_seqs[header]))
        # print(sum([a,t,c,g,n]))

        # replacements = [(b, "") for b in "ATCGNatcgn"];
        # tmp = cur_seqs[header];
        # for pat, repl in replacements:
        #     tmp = re.sub(pat, repl, tmp);

        # print(header);
        # print(tmp);
        # print("----")

        assert seqlen == sum([a,t,c,g,n]), " * ERROR: Non-standard bases present: " + consensus_file + " " + header;

        perc_n = n / seqlen;

        if n > 5000:
            n_filter = "FILTER";
        else:
            n_filter = "PASS";

        outline = [ str(val) for val in [sample, header, a, t, c, g, n, seqlen, perc_n, n_filter] ];
        outfile.write(",".join(outline) + "\n");

############################################################