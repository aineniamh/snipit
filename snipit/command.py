from snipit import __version__
import argparse
import os
import pkg_resources
from . import _program
from Bio import SeqIO
import collections
import snp_functions as sfunks
import sys

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='snipit', 
    usage='''snipit <alignment> [options]''')

    parser.add_argument('alignment',help="Input alignment fasta file")
    parser.add_argument('-o',"--output-dir",action="store",dest="output_dir")
    parser.add_argument("-r","--reference", action="store",help="Indicates which sequence in the alignment is the reference", dest="reference")
    parser.add_argument("-f","--format",action="store",default="png")
    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    lengths = []
    lengths_info = []
    num_seqs = 0
    for record in SeqIO.parse(args.alignment, "fasta"):
        lengths.append(len(record))
        lengths_info.append((record.id, len(record)))
        num_seqs +=1
    print(f"{num_seqs} found in alignment file")
    if len(set(lengths))!= 1:
        sys.stderr.write("Error: not all of the sequences in the alignment are the same length\n")
        for i in lengths_info:
            print(f"{i[0]}\t{i[1]}\n")
        sys.exit(-1)
    
    if not args.output_dir:
        output_dir = cwd
    else:
        output_dir = os.path.join(cwd, args.output_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    snp_file = os.path.join(output_dir,"snps.csv")
    ambiguities_file = os.path.join(output_dir,"ambiguities.csv")
    output = os.path.join(output_dir,f"genome_graph.{args.format}")
    
    reference,alignment = sfunks.get_ref_and_alignment(args.alignment,args.reference)

    snp_dict,record_snps = sfunks.find_snps(reference,alignment)

    record_ambs = sfunks.find_ambiguities(alignment, snp_dict)

    sfunks.make_graph(num_seqs,snp_file,record_ambs,record_snps,output)

if __name__ == '__main__':
    main()

    