import sys
import argparse

#parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="extract longest contig from a fasta file")
    parser.add_argument("-i", "--input",
    help="input fasta file (contigs.fasta)",
    required=True)
    parser.add_argument("-o", "--output",
    help="output fasta file (longest contig)",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#keep track of the current contig header and sequence
current_header = None
current_seq_parts = []

#keep track of the best(longest) contig so far
best_header = None
best_seq = ""
best_len = 0

#open and read fasta file
with open(infile, "r") as f:
    for line in f:
        line = line.strip()

        #skip empty lines
        if line == "":
            continue

        #header line starts a new contig
        if line.startswith(">"):
            #if we were already building a contig, finalize it
            if current_header is not None:
                #join sequence lines into one string
                seq = "".join(current_seq_parts)
                seq_len = len(seq)

                #update best contig if this one is longer
                if seq_len > best_len:
                    best_len = seq_len
                    best_header = current_header
                    best_seq = seq

            #start new contig
            current_header = line
            current_seq_parts = []

        else:
            #sequence line, append to list
            current_seq_parts.append(line)

#after loop ends, finalize the last contig
if current_header is not None:
    seq = "".join(current_seq_parts)
    seq_len = len(seq)
    if seq_len > best_len:
        best_len = seq_len
        best_header = current_header
        best_seq = seq

#write the best contig to output fasta
with open(outfile, "w") as out:
    #write header
    out.write(best_header + "\n")

    #write sequence wrapped to 60 chars per line for readability
    width = 60
    for i in range(0, len(best_seq), width):
        out.write(best_seq[i:i+width] + "\n")