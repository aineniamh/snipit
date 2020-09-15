#!/usr/bin/env python3
import os
import argparse
import collections
from Bio import SeqIO
import csv

cwd = os.getcwd()

def parse_args():
    parser = argparse.ArgumentParser(description='Find ambiguities at snp sites.')

    parser.add_argument("--input", action="store", type=str, dest="input")
    parser.add_argument("--report", action="store", type=str, dest="report")
    parser.add_argument("--output", action="store", type=str, dest="output")

    return parser.parse_args()

def find_snps(input_file, output, report):

    input_seqs = collections.defaultdict(list)
    outgroup_seq = ""

    for record in SeqIO.parse(input_file, "fasta"):
        if record.id == "outgroup":
            outgroup_seq = record.seq.upper()
        else:
            input_seqs[str(record.seq).upper()].append(record.id)

    with open(output, "w") as fw:
        fw.write("name\tnum_snps\tambiguous_snps\n")
        snp_dict = collections.defaultdict(list)
        with open(report, newline="") as f:
            reader = csv.DictReader(f, delimiter="\t")
            for row in reader:
                snps = row["snps"].split(";")
                for snp in snps:
                    x_position = int(snp[:-2])-1
                    ref = snp[-2]
                    snp_dict[x_position]=ref

        for query_seq in input_seqs:
            snps =[]

            for i in snp_dict:
                bases = [query_seq[i],outgroup_seq[i]]
                if bases[0] != bases[1]:
                    if bases[0] not in ["A","T","G","C"]:
                        
                        snp = f"{i+1}{bases[1]}{bases[0]}" # position-outgroup-query
                        snps.append(snp)
            
            for record_name in input_seqs[query_seq]:
                print(record_name, len(snps))
                snp_str = ";".join(snps)
                fw.write(f"{record_name}\t{len(snps)}\t{snp_str}\n")

if __name__ == '__main__':

    args = parse_args()

    find_snps(args.input, args.output, args.report)