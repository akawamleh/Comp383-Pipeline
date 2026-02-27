import sys
import argparse
import gzip

#parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="count read pairs from R1 fastq or fastq.gz")
    parser.add_argument("-i","--input",help="input R1 fastq(.gz)",required=True)
    parser.add_argument("-o","--output",help="output file",required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#choose open function based on extension
def open_maybe_gzip(path):
    if path.endswith(".gz"):
        return gzip.open(path,"rt")
    return open(path,"r")

#line_count/4=number_of_reads
line_count = 0
with open_maybe_gzip(infile) as f:
    for _ in f:
        line_count += 1

pairs = line_count // 4

with open(outfile,"w") as out:
    out.write(str(pairs) + "\n")