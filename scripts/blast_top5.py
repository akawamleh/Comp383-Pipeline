import sys
import argparse

#parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="keep top 5 blast hits and enforce one row per subject accession")
    parser.add_argument("-i","--input",help="input blast tsv (outfmt 6)",required=True)
    parser.add_argument("-o","--output",help="output tsv (top 5)",required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
infile = arguments.input
outfile = arguments.output

#store best row per subject accession
best_by_sacc = {}

#read blast output
with open(infile,"r") as f:
    for line in f:
        line = line.rstrip("\n")
        if line == "":
            continue

        parts = line.split("\t")

        #expected columns:
        #0sacc 1pident 2length 3qstart 4qend 5sstart 6send 7bitscore 8evalue 9stitle
        sacc = parts[0]
        bitscore = float(parts[7])
        evalue = float(parts[8])

        #keep the best row for each sacc using bitscore, break ties with evalue
        if sacc not in best_by_sacc:
            best_by_sacc[sacc] = (bitscore, evalue, parts)
        else:
            old_bitscore, old_evalue, _ = best_by_sacc[sacc]
            if (bitscore > old_bitscore) or (bitscore == old_bitscore and evalue < old_evalue):
                best_by_sacc[sacc] = (bitscore, evalue, parts)

#sort all unique subject hits by bitscore desc, then evalue asc
rows = list(best_by_sacc.values())
rows.sort(key=lambda x: (-x[0], x[1]))

#take top 5
top5 = rows[:5]

#write output
with open(outfile,"w") as out:
    for bitscore, evalue, parts in top5:
        out.write("\t".join(parts) + "\n")