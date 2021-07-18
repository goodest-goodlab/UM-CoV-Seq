#!/usr/bin/python
############################################################
# Replaces bar code ids with sample ids in consensus files
############################################################

import sys, os, csv, datetime

############################################################
# Functions

def getDateTime():
#Function to get the date and time in a certain format.
	return datetime.datetime.now().strftime("%m.%d.%Y | %H:%M:%S");

###############

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

input_dir, batch, sample_file, summary_file, output_file = sys.argv[1:];
# Parse positional input arguments:
# input_dir: The directory containing individual consensus sequences for each sample
# batch: The batch/run ID
# sample_file: The main sample table with barcodes and sample IDs
# summary_file: The output table from compile_results.r
# output_file: Path to write the output sequences.

outdir = os.path.dirname(output_file);
print(outdir);
if not os.path.isdir(outdir):
    os.system("mkdir " + outdir);
# Create the gisaid directory if it doesn't exist... not sure why snakemake wasn't doing this

print("# " + getDateTime() + " Reading sample information: " + sample_file);
samples = {};
# samples will be the key between the barcode and sample ID, as samples[<barcode>] = <sample id>

icv_list = ["F3189"];
# Keep track of that one sample that ends in "-ICV"...

first = True;
for line in open(sample_file):
# Read the sample file
    line = line.strip().split(",");
    # Split the lines by commas
    if first:
        meta_headers = line;

        barcode_ind = line.index("Barcode");
        sample_ind = line.index("ID By UMGC");
        batch_ind = line.index("MiSeq Run");
        first = False;
        continue;
    # For the first line, get the indices of the barcode, sample ID, and batch ID

    if line[batch_ind] == batch:
        samples[line[barcode_ind]] = line[sample_ind];
    # Save the barcode and sample id if the current sample is in the current batch

print("# " + getDateTime() + " " + str(len(samples)) + " samples read");

###############

print("# " + getDateTime() + " Reading summary information: " + summary_file);
with open(summary_file) as csv_file:
    csv_reader = csv.reader(csv_file, delimiter=',', quotechar='"');
    # Read the summary file as a csv file

    summary = {};
    # The summary info retrieved from the main results table will simply be whether the lineage predictions match
    first = True;
    for row in csv_reader:
        if first:
            barcode_ind = row.index("sample");
            pango_match_ind = row.index("pango.match");
            nextclade_match_ind = row.index("nextclade.match");
            first = False;
            continue;
        # For the first line in the summary file, get the indices of the barcode and lineage matches
        # Lineage matches are encoded as "1" for yes and "0" for no

        summary[row[barcode_ind]] = [ row[pango_match_ind], row[nextclade_match_ind] ];
        # Save info for current sample as list in form of [ <pangolin match>, <nextclade match> ]

print("# " + getDateTime() + " " + str(len(samples)) + " samples read");

###############

print("# " + getDateTime() + " Relabeling samples and writing to file: " + output_file);
with open(output_file, "w") as outfile:
    for barcode in samples:
        if barcode not in summary:
            print(" * WARNING: Cannot find sample in summary file. Skipping: " + barcode);
            continue;

        if summary[barcode] != ["1", "1"]:
            continue;
        # For every sample, check if the lineage predictions match. If not, do not write thie sample
        # to the GISAID fasta file

        sample = samples[barcode];
        # Get the sample ID from the barcode
        
        if barcode in icv_list:
            barcode += "-ICV";
        # For that one sample...

        fafile = os.path.join(input_dir, barcode + ".fa");
        assert os.path.isfile(fafile), "* ERROR: Cannot find fasta file: " + fafile;
        # Get the filename of the sample's consensus sequence and exit if it is not found

        cur_seq, cur_titles = readFasta(fafile);
        # Read the current sample consensus sequence

        gisaid_title = ">hCoV-19/USA/MT-" + sample + "/2021";
        # The new header for GISAID using the UMGC sample ID

        outfile.write(gisaid_title + "\n");
        outfile.write(cur_seq[cur_titles[0]] + "\n");
        # Write the current sequence

############################################################

