#!/usr/bin/env python3

# imports of built-ins
import sys
import os
import argparse
import pkg_resources
import collections

# imports from other modules
from Bio import SeqIO

# imports from this module
from snipit import __version__
from . import _program
from snipit.scripts import snp_functions as sfunks

thisdir = os.path.abspath(os.path.dirname(__file__))
cwd = os.getcwd()

def bp_range(s):
    """
        Crude function to parse positions or position ranges (inclusive) passed as a string by argparse.
        Input: string in the format "100-200" or "100"
        Returns a list with integer positions.
        Arguably better solved by a regex, but still would need to typecast
    """
    # try to parse as a range
    try:
        start,end = map(int, s.split('-'))
        return list(range(start,end+1))
    except ValueError:
        # if range parsing fails, perhaps it's only one position. try to parse as a single int
        try: 
            pos = int(s)
            return [pos]
        except ValueError:
            raise argparse.ArgumentTypeError("Coordinates must be in the format 'start-end' or 'pos'")
        

def main(sysargs = sys.argv[1:]):

    parser = argparse.ArgumentParser(prog = _program, 
    description='snipit', 
    usage='''snipit <alignment> [options]''')

    parser.add_argument('alignment',help="Input alignment fasta file")
    parser.add_argument("-r","--reference", action="store",help="Indicates which sequence in the alignment is\nthe reference (by sequence ID).\nDefault: first sequence in alignment", dest="reference")
    parser.add_argument("-l","--labels", action="store",help="Optional csv file of labels to show in output snipit plot. Default: sequence names", dest="labels")
    parser.add_argument("--l-header", action="store",help="Comma separated string of column headers in label csv. First field indicates sequence name column, second the label column. Default: 'name,label'", dest="label_headers",default="name,label")


    parser.add_argument('-d',"--output-dir",action="store",help="Output directory. Default: current working directory", dest="output_dir")
    parser.add_argument('-o',"--output-file",action="store",help="Output file name stem. Default: snp_plot", default="snp_plot",dest="outfile")
    parser.add_argument('-s',"--write-snps",action="store_true",help="Write out the SNPs in a csv file.",dest="write_snps")
    parser.add_argument("-f","--format",action="store",help="Format options (png, jpg, pdf, svg, tiff) Default: png",default="png")

    parser.add_argument("--height",action="store",type=float,help="Overwrite the default figure height",default=0)
    parser.add_argument("--width",action="store",type=float,help="Overwrite the default figure width",default=0)
    parser.add_argument("--size-option",action="store",help="Specify options for sizing. Options: expand, scale",dest="size_option",default="scale")
    parser.add_argument("--solid-background",action="store_true",help="Force the plot to have a solid background, rather than a transparent one.",dest="solid_background")
    parser.add_argument("--flip-vertical",action='store_true',help="Flip the orientation of the plot so sequences are below the reference rather than above it.",dest="flip_vertical")

    parser.add_argument('--include-positions', dest='included_positions', type=bp_range, nargs='+', default=None, help="One or more range (closed, inclusive; one-indexed) or specific position only included in the output. Ex. '100-150' or Ex. '100 101' Considered before '--exclude-positions'.")
    parser.add_argument('--exclude-positions', dest='excluded_positions', type=bp_range, nargs='+', default=None, help="One or more range (closed, inclusive; one-indexed) or specific position to exclude in the output. Ex. '100-150' or Ex. '100 101' Considered after '--include-positions'.")
    parser.add_argument("--exclude-ambig-pos",dest="exclude_ambig_pos",action='store_true',help="Exclude positions with ambig base in any sequences. Considered after '--include-positions'")
    parser.add_argument("--sort-by-mutation-number", action='store_true',
                        help="Render the graph with sequences sorted by the number of SNPs relative to the reference (fewest to most). Default: False", dest="sort_by_mutation_number")
    parser.add_argument("--sort-by-id", action='store_true',
                        help="Sort sequences alphabetically by sequence id. Default: False", dest="sort_by_id")
    parser.add_argument("--sort-by-mutations", type=str, help="Sort sequences by bases at specified positions. Positions are comma separated integers. Ex. '1,2,3'", dest="sort_by_mutations")
    parser.add_argument("--high-to-low", action='store_false',
                        help="If sorted by mutation number is selected, show the sequences with the fewest SNPs closest to the reference. Default: False",
                        dest="high_to_low")

    parser.add_argument("-c","--colour-palette",dest="colour_palette",action="store",help="Specify colour palette. Options: primary, classic, purine-pyrimidine, greyscale, wes, verity",default="classic")
    parser.add_argument("--recombi-mode",action='store_true',dest="recombi_mode",help="Allow colouring of query seqeunces by mutations present in two 'recombi-references' from the input alignment fasta file")
    parser.add_argument("--recombi-references",action='store',type=str,dest="recombi_references",help="Specify two comma separated sequence IDs in the input alignment to use as 'recombi-references'. Ex. Sequence_ID_A,Sequence_ID_B")

    if len(sysargs)<1:
        parser.print_help()
        sys.exit(-1)
    else:
        args = parser.parse_args(sysargs)

    num_seqs,ref_input,record_ids,length = sfunks.qc_alignment(args.alignment,args.reference,cwd)
        
    
    if args.reference:
        ref_file,ref_input = sfunks.reference_qc(args.reference, record_ids,cwd)
    else:
        sfunks.check_ref(args.recombi_mode)

    if args.recombi_references:
        sfunks.recombi_qc(args.recombi_references, args.reference, record_ids,cwd)

    if not args.output_dir:
        output_dir = cwd
    else:
        output_dir = os.path.join(cwd, args.output_dir)
        if not os.path.exists(output_dir):
            os.mkdir(output_dir)

    output = os.path.join(output_dir,f"{args.outfile}.{args.format}")
    
    label_map = sfunks.label_map(record_ids,args.labels,args.label_headers,cwd)

    reference,alignment = sfunks.get_ref_and_alignment(args.alignment,ref_input,label_map)

    snp_dict,record_snps,num_snps = sfunks.find_snps(reference,alignment)

    record_ambs = sfunks.find_ambiguities(alignment, snp_dict)

    colours = sfunks.get_colours(args.colour_palette)

    sfunks.check_format(args.format)
    sfunks.check_size_option(args.size_option)

    sfunks.write_out_snps(args.write_snps,record_snps,output_dir)

    sfunks.make_graph(num_seqs,
                        num_snps,
                        record_ambs,
                        record_snps,
                        output,
                        label_map,
                        colours,
                        length,
                        args.width,
                        args.height,
                        args.size_option,
                        args.solid_background,
                        args.flip_vertical,
                        args.included_positions,
                        args.excluded_positions,
                        args.exclude_ambig_pos,
                      args.sort_by_mutation_number,
                      args.high_to_low,
                      args.sort_by_id,
                      args.sort_by_mutations,
                      args.recombi_mode,
                      args.recombi_references)


if __name__ == '__main__':
    main()
