import sys
import argparse

#parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="compute contigs >1000 bp and total bp")
    parser.add_argument("-i", "--input",
    help="input fasta file",
    required=True)
    parser.add_argument("-o", "--output",
    help="output file",
    required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#threshold for contig length
MIN_LEN = 1000

#initialize counters
num_over_1000 = 0
total_bp_over_1000 = 0
current_length = 0

#read fasta file
with open(infile, "r") as f:
    for line in f:
        line = line.strip()

        #skip empty lines
        if line == "":
            continue

        #new contig header
        if line.startswith(">"):
            #check previous contig
            if current_length > 0:
                if current_length > MIN_LEN:
                    num_over_1000 += 1
                    total_bp_over_1000 += current_length
            #reset length
            current_length = 0
        else:
            #add sequence length
            current_length += len(line)

#check last contig
if current_length > 0:
    if current_length > MIN_LEN:
        num_over_1000 += 1
        total_bp_over_1000 += current_length

#write output
with open(outfile, "w") as out:
    out.write(str(num_over_1000) + " " + str(total_bp_over_1000) + "\n")