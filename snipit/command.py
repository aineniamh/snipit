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
    parser.add_argument("-r","--reference", action="store",help="Indicates which sequence in the alignment is\nthe reference (by sequence ID).\nCan also accept a genbank file.\nDefault: first sequence in alignment", dest="reference")

    parser.add_argument('-o',"--output-dir",action="store",help="Output directory. Default: current working directory", dest="output_dir")
    parser.add_argument("-f","--format",action="store",help="Format options (png, jpg, pdf, svg, tiff). Default: png",default="png")

    parser.add_argument("--height",action="store",type=float,help="Overwrite the default figure height",default=0)
    parser.add_argument("--width",action="store",type=float,help="Overwrite the default figure width",default=0)
    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    num_seqs,ref_input,record_ids,length = sfunks.qc_alignment(args.alignment,args.reference,cwd)
        
    
    if args.reference:
        ref_file,ref_input = sfunks.reference_qc(args.reference, record_ids,cwd)

    if not args.output_dir:
        output_dir = cwd
    else:
        output_dir = os.path.join(cwd, args.output_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    output = os.path.join(output_dir,f"genome_graph.{args.format}")
    
    reference,alignment = sfunks.get_ref_and_alignment(args.alignment,ref_input)

    snp_dict,record_snps,num_snps = sfunks.find_snps(reference,alignment)

    record_ambs = sfunks.find_ambiguities(alignment, snp_dict)

    sfunks.make_graph(num_seqs,num_snps,record_ambs,record_snps,output,length,args.width,args.height)

if __name__ == '__main__':
    main()

    