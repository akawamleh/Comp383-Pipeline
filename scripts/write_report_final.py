import sys
import argparse

#parse command line arguments
def check_arg(args=None):
    parser = argparse.ArgumentParser(description="write PipelineReport.txt from pipeline outputs")
    parser.add_argument("-o","--output",help="PipelineReport.txt",required=True)
    parser.add_argument("-s","--samples",help="sample names",nargs="+",required=True)
    parser.add_argument("-r","--results",help="results directory for this mode",required=True)
    return parser.parse_args(args)

#retrieve command line arguments
arguments = check_arg(sys.argv[1:])
outfile = arguments.output
samples = arguments.samples
results_dir = arguments.results.rstrip("/")

with open(outfile,"w") as out:
    for s in samples:
        before_path = f"{results_dir}/counts/{s}_before.txt"
        after_path = f"{results_dir}/counts/{s}_after.txt"
        stats_path = f"{results_dir}/assemblies/{s}/contigs_over_1000_stats.txt"
        top5_path = f"{results_dir}/assemblies/{s}/{s}_blast_top5.tsv"

        with open(before_path,"r") as f:
            before_pairs = f.read().strip()

        with open(after_path,"r") as f:
            after_pairs = f.read().strip()

        with open(stats_path,"r") as f:
            n_over, total_over = f.read().strip().split()

        out.write(f"Before Bowtie2 filtering, {s} has {before_pairs} read pairs.\n")
        out.write(f"After Bowtie2 filtering, {s} has {after_pairs} read pairs.\n")
        out.write(f"In the assembly of sample {s}, there are {n_over} contigs > 1000 bp and {total_over} total bp.\n")
        out.write(f"Sample {s} top 5 blast hits:\n")
        out.write("sacc\tpident\tlength\tqstart\tqend\tsstart\tsend\tbitscore\tevalue\tstitle\n")

        with open(top5_path,"r") as f:
            for line in f:
                out.write(line if line.endswith("\n") else line + "\n")

        out.write("\n")