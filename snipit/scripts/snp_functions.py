#!/usr/bin/env python3

# imports of built-ins
import os
import sys
import argparse
import collections
from itertools import cycle, chain
import csv
import math
from itertools import groupby, count
from collections import OrderedDict

# imports from other modules
from Bio import SeqIO
from Bio.Seq import Seq
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.patches as patches
from matplotlib.patches import Polygon

colour_list = ["lightgrey","white"]
colour_cycle = cycle(colour_list)
END_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'


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
        


def check_ref(recombi_mode):
    if recombi_mode:
        sys.stderr.write(red(f"Error: Please explicitly state reference sequence when using `--recombi-mode`\n"))
        sys.exit(-1)


def qc_alignment(alignment,reference,cds_mode,cwd):
    lengths = []
    lengths_info = []
    num_seqs = 0

    record_ids = []
    ref_input = ""

    alignment_file = os.path.join(cwd, alignment)
    if not os.path.exists(alignment_file):
        sys.stderr.write(red(f"Error: can't find alignment file at {alignment_file}\n"))
        sys.exit(-1)

    try:
        for record in SeqIO.parse(alignment_file, "fasta"):
            if ref_input == "":
                ref_input = record.id
            lengths.append(len(record))
            record_ids.append(record.id)
            lengths_info.append((record.id, len(record)))
            num_seqs +=1
    except:
        sys.stderr.write(red(f"Error: alignment file must be in fasta format\n"))
        sys.exit(-1)

    if num_seqs == 1:
        if reference:
            if reference.split(".")[-1] not in ["gb","genbank"]:
                sys.stderr.write(red(f"Error: alignment file must contain more than just the reference. Either provide a reference genbank file or add more sequences to your alignment.\n"))
                sys.exit(-1)
        else:
            sys.stderr.write(red(f"Error: alignment file must contain more than just the reference. Either provide a reference genbank file or add more sequences to your alignment.\n"))
            sys.exit(-1)
    unique_lengths = set(lengths)
    if len(unique_lengths)!= 1:
        sys.stderr.write(red("Error: not all of the sequences in the alignment are the same length\n"))
        for i in lengths_info:
            print(f"{i[0]}\t{i[1]}\n")
        sys.exit(-1)
    
    if cds_mode and unique_lengths[0]%3!=0:
        sys.stderr.write(red("Error: CDS mode flag used but alignment length not a multiple of 3.\n"))
        sys.exit(-1)

    return num_seqs,ref_input,record_ids,lengths[0]

def reference_qc(reference, record_ids,cwd):
    ref_file = ""
    if "." in reference and reference.split(".")[-1] in ["gb","genbank"]:
        ref_path = os.path.join(cwd, reference)
        if not os.path.exists(ref_path):
            sys.stderr.write(red(f"Error: can't find genbank file at {ref_path}\n"))
            sys.exit(-1)
        else:
            ref_input = ""
            ref_file = ""
            record_count = 0
            for record in SeqIO.parse(reference, "genbank"):
                ref_input = record.seq
                ref_file = record
                record_count += 1

            if record_count >1:
                sys.stderr.write(red(f"Error: more than one record found in reference genbank file\n"))
                sys.exit(-1)

    elif reference not in record_ids:
        sys.stderr.write(red(f"Error: input reference {reference} not found in alignment\n"))
        sys.exit(-1)

    else:
        ref_input = reference

    return ref_file, ref_input


def recombi_qc(recombi_refs, reference, record_ids,cwd):
    recombi_refs = recombi_refs.split(",")
    if not len(recombi_refs) == 2:
        sys.stderr.write(red(f"Error: input 2 references separated by a comma for `--recombi-references`.\n"))
        sys.exit(-1)

    for ref in recombi_refs:
        if ref == "":
            sys.stderr.write(red(f"Error: input 2 references separated by a comma for `--recombi-references`.\n"))
            sys.exit(-1)
        if ref == reference:
            sys.stderr.write(red(f"Error: please input a distinct outgroup reference from the parent recombinant references specified in `--recombi-references`.\n"))
            sys.exit(-1)
        if ref not in record_ids:
            sys.stderr.write(red(f"Error: please check references specified in `--recombi-references` match a sequence name in the input alignment.\n"))
            sys.exit(-1)



def label_map(record_ids,labels,column_names,cwd):
    seq_col,label_col = column_names.split(",")

    label_map = {}
    if labels:
        label_file = os.path.join(cwd,labels)
        if os.path.exists(label_file):
            with open(label_file, "r") as f:
                reader = csv.DictReader(f)
                if seq_col not in reader.fieldnames:
                    sys.stderr.write(red(f"Error: {seq_col} not a column name in {labels}\n"))
                    sys.exit(-1)
                elif label_col not in reader.fieldnames:
                    sys.stderr.write(red(f"Error: {label_col} not a column name in {labels}\n"))
                    sys.exit(-1)

                for row in reader:
                    sequence_name,label = (row[seq_col],row[label_col])
                    label_map[sequence_name] = label
        else:
            sys.stderr.write(red(f"Error: {record_id} not in {labels} header\n"))
            sys.exit(-1)

        for record_id in record_ids:
            if record_id not in label_map:
                sys.stderr.write(red(f"Error: {record_id} not in {labels} header\n"))
                sys.exit(-1)
    else:
        for record_id in record_ids:
            label_map[record_id] = record_id

    return label_map

def next_colour():
    return next(colour_cycle)

def get_ref_and_alignment(input_file,reference,label_map):
    input_seqs = collections.defaultdict(list)
    reference_seq = ""

    for record in SeqIO.parse(input_file, "fasta"):
        if record.id == reference:
            reference_seq = record.seq.upper()
            if record.id not in label_map:
                label_map["reference"]=record.id
            else:
                label_map["reference"]=label_map[record.id]
        else:
            input_seqs[str(record.seq).upper()].append(record.id)

    return reference_seq, input_seqs

def merge_indels(indel_list,prefix):
    if indel_list:
        groups = groupby(indel_list, key=lambda item, c=count():item-next(c))
        tmp = [list(g) for k, g in groups]
        merged_indels = []
        for i in tmp:
            indel = f"{i[0]}:{prefix}{len(i)}"
            merged_indels.append(indel)
        return merged_indels

    return indel_list

def find_snps(reference_seq,input_seqs,show_indels):
    non_amb = ["A","T","G","C"]
    snp_dict = {}

    record_snps = {}
    var_counter = collections.Counter()
    for query_seq in input_seqs:
        snps =[]
        insertions = []
        deletions = []
        for i in range(len(query_seq)):
            bases = [query_seq[i],reference_seq[i]]
            if bases[0] != bases[1]:
                if bases[0] in non_amb and bases[1] in non_amb:

                    snp = f"{i+1}:{bases[1]}{bases[0]}" # position-reference-query

                    snps.append(snp)
                elif bases[0]=='-' and show_indels:
                #if there's a gap in the query, means a deletion
                    deletions.append(i+1)
                elif bases[1]=='-' and show_indels:
                    #if there's a gap in the ref, means an insertion
                    insertions.append(i+1)

        if show_indels:
            insertions = merge_indels(insertions,"ins")
            deletions = merge_indels(deletions,"del")

        variants = []
        for var_list in [snps,insertions,deletions]:
            for var in var_list:
                var_counter[var]+=1
                variants.append(var)

        variants = sorted(variants, key = lambda x : int(x.split(":")[0]))

        snp_dict[query_seq] = variants

        for record in input_seqs[query_seq]:
            record_snps[record] = variants

    return snp_dict,record_snps,len(var_counter)

def find_ambiguities(alignment, snp_dict):

    snp_sites = collections.defaultdict(list)
    for seq in snp_dict:
        snps = snp_dict[seq]
        for snp in snps:
            pos,var = snp.split(":")
            index = int(pos)-1

            ref_allele = var[0]
            snp_sites[index]=ref_allele

    amb_dict = {}

    for query_seq in alignment:
        snps =[]

        for i in snp_sites:
            bases = [query_seq[i],snp_sites[i]] #if query not same as ref allele
            if bases[0] != bases[1]:
                if bases[0] not in ["A","T","G","C"]:

                    snp = f"{i+1}:{bases[1]}{bases[0]}" # position-outgroup-query
                    snps.append(snp)

        for record in alignment[query_seq]:
            amb_dict[record] = snps

    return amb_dict


def recombi_ref_snps(recombi_references, snp_records):

    #print(recombi_references)
    recombi_refs = recombi_references.split(",")
    recombi_snps = []
    #print(recombi_refs)
    for ref in recombi_refs:
        recombi_snps.append(snp_records[ref])

    #print(recombi_snps)
    return recombi_snps,recombi_refs

def recombi_painter(snp_to_check,recombi_snps):

    recombi_ref_1 = recombi_snps[0]
    recombi_ref_2 = recombi_snps[1]
    common_snps = []

    for snp in recombi_ref_1:
        if snp in recombi_ref_2:
            common_snps.append(snp)

    if snp_to_check in common_snps:
        return "Both"
    elif snp_to_check in recombi_ref_1:
        return "lineage_1"
    elif snp_to_check in recombi_ref_2:
        return "lineage_2"
    else:
        return "Private"


def write_out_snps(write_snps,record_snps,output_dir):
    with open(os.path.join(output_dir,"snps.csv"),"w") as fw:
        fw.write("record,snps,num_snps\n")
        for record in record_snps:
            snps = ";".join(record_snps[record])
            fw.write(f"{record},{snps},{len(record_snps[record])}\n")


def make_graph(num_seqs,
                num_snps,
                amb_dict,
                snp_records,
                output,
                label_map,
                colour_dict,
                length,
                width,
                height,
                size_option,
                solid_background,
               flip_vertical=False,
               included_positions=None,
               excluded_positions=None,
               exclude_ambig_pos=False,
               sort_by_mutation_number=False,
               high_to_low=True,
               sort_by_id=False,
               sort_by_mutations=False,
               recombi_mode=False,
               recombi_references=[]
               ):
    y_level = 0
    ref_vars = {}
    snp_dict = collections.defaultdict(list)
    included_positions = set(chain.from_iterable(included_positions)) if included_positions is not None else set()
    excluded_positions = set(chain.from_iterable(excluded_positions)) if excluded_positions is not None else set()

    if sort_by_mutation_number:
        snp_counts = {}
        for record in snp_records:
            snp_counts[record] = int(len(snp_records[record]))
        ordered_dict = dict(sorted(snp_counts.items(), key=lambda item: item[1], reverse=high_to_low))
        record_order = list(OrderedDict(ordered_dict).keys())

    elif sort_by_id:
        record_order = list(sorted(snp_records.keys()))

    elif sort_by_mutations:
        mutations = sort_by_mutations.split(",")
        sortable_record = {}
        for record in snp_records:
            bases = []
            for sort_mutation in mutations:
                found = False
                for record_mutation in snp_records[record]:
                    if int(record_mutation.split(":")[0]) == int(sort_mutation):
                        bases.append(record_mutation[-1])
                        found = True
                        break
                if not found:
                    bases.append("0")
            sortable_record[record] = "".join(bases) + record
        record_order = list(OrderedDict(sorted(sortable_record.items(), key=lambda item: item[1], reverse=high_to_low)).keys())

    else:
        record_order = list(snp_records.keys())

    if recombi_mode:
        # Get a list of SNPs present in each recombi_reference
        recombi_snps,recombi_refs = recombi_ref_snps(recombi_references, snp_records)
        # Set the colour palette to "recombi"
        colour_dict = get_colours("recombi")
        # Reorder list to put recombi_references at the start
        record_order.remove(recombi_refs[0])
        record_order.insert(0, recombi_refs[0])
        record_order.remove(recombi_refs[1])
        record_order.insert(1, recombi_refs[1])

    for record in record_order:

        # y level increments per record, add a gap after the two recombi_refs
        if recombi_mode and y_level == 2:
            y_level += 1.2
        else:
            y_level +=1

        # for each record get the list of snps
        snps = snp_records[record]
        x = []
        y = []

        for snp in snps:
            # snp => 12345AT
            pos,var = snp.split(":")
            x_position = int(pos)
            if var.startswith("del"):
                length_indel = var[3:]
                ref = f"{length_indel}"
                base = "-"
            elif var.startswith("ins"):
                length_indel = var[3:]
                ref = "-"
                base = f"{length_indel}"
            else:
                ref = var[0]
                base = var[1]

            ref_vars[x_position]=ref
            if recombi_mode:
                recombi_out = recombi_painter(snp, recombi_snps)
                # Add name of record, ref, SNP in record, y_level, if SNP is in either recombi_reference...
                snp_dict[x_position].append((record, ref, base, y_level, recombi_out))
            else:
                # ...otherwise add False instead to help the colour logic
                snp_dict[x_position].append((record, ref, base, y_level, False))

        # if there are ambiguities in that record, add them to the snp dict too
        if record in amb_dict:
            for amb in sorted(amb_dict[record]):
                # amb => 12345AN
                pos,var = amb.split(":")
                x_position = int(pos)

                # if positions with any ambiguities should be ignored, note the position
                if exclude_ambig_pos:
                    excluded_positions.add(x_position)
                else:
                    ref = var[0]
                    base = var[1]
                    ref_vars[x_position]=ref
                    # Add name of record, ref, SNP in record, y_level and False for "recombi_mode" colour logic
                    snp_dict[x_position].append((record, ref, base, y_level, False))

    # gather the positions that are not explicitly excluded,
    # but are not among those to be included
    positions_not_included=set()
    if len(included_positions)>0:
        # of the positions present,
        # gather a set of positions which should NOT be included in the output
        positions_not_included = set(snp_dict.keys()) - included_positions

    # remove positions which should be ignored or are not included (pop items from union of the two sets)
    for pos in excluded_positions | positions_not_included:
        # remove records for the position, if present
        snp_dict.pop(pos, None)

    spacing = length/(len(snp_dict)+1)
    y_inc = (spacing*0.8*y_level)/length

    if size_option == "expand":
        if not width:
            if num_snps ==0:
                print(red(f"Note: no SNPs found between the reference and the alignment"))
                width = 10
            else:
                if len(snp_dict) <10:
                    width = 10
                else:
                    width = 0.25* len(snp_dict)

        if not height:
            if y_level < 5:
                height = 5
            else:
                height = (y_inc*3 + 0.5*y_level + y_inc*2) # bottom chunk, and num seqs, and text on top

    elif size_option == "scale":
        if not width:
            if num_snps == 0:
                print(red(f"Note: no SNPs found between the reference and the alignment"))
                width = 12
            else:
                width = math.sqrt(num_snps)*3

        if not height:
            height = math.sqrt(num_seqs)*2
            y_inc = 1

    # if the plot is flipped vertically, place the x-axis (genome map) labels on top
    if flip_vertical:
        plt.rcParams['xtick.bottom'] = plt.rcParams['xtick.labelbottom'] = False
        plt.rcParams['xtick.top'] = plt.rcParams['xtick.labeltop'] = True

    # width and height of the figure
    fig, ax = plt.subplots(1,1, figsize=(width,height), dpi=250)

    y_level = 0

    for record in record_order:

        # y position increments, with a gap after the two recombi_refs
        if recombi_mode and y_level == 2:
            y_level += y_inc + 0.2
        else:
            y_level += y_inc


        # either grey or white
        col = next_colour()

        # for each record (sequence) draw a rectangle the length of the whole genome (either grey or white)
        rect = patches.Rectangle((0,y_level-(0.5*y_inc)), length, y_inc ,alpha=0.3, fill=True, edgecolor='none',facecolor=col)
        ax.add_patch(rect)

        # for each record add the name to the left hand side
        ax.text(-50, y_level, label_map[record], size=9, ha="right", va="center")

    position = 0

    for snp in sorted(snp_dict):
        position += spacing

        # write text adjacent to the SNPs shown with the numeric position
        # the text alignment is toggled right/left (top/bottom considering 90-deg rotation) if the plot is flipped
        ax.text(position, y_level+(0.55*y_inc), snp, size=9, ha="center", va="bottom" if not flip_vertical else "top", rotation=90)

        # snp position labels
        left_of_box = position-(0.4*spacing)
        right_of_box = position+(0.4*spacing)

        top_polygon = y_inc * -0.7
        bottom_polygon = y_inc * -1.7

        for sequence in snp_dict[snp]:

            name,ref,var,y_pos,recombi_out = sequence
            bottom_of_box = (y_pos*y_inc)-(0.5*y_inc)
            # draw box for snp
            if recombi_out:
                rect = patches.Rectangle((left_of_box,bottom_of_box),spacing*0.8,  y_inc,alpha=0.5, fill=True, edgecolor='none',facecolor=colour_dict[recombi_out])
            elif var in colour_dict:
                rect = patches.Rectangle((left_of_box,bottom_of_box),spacing*0.8,  y_inc,alpha=0.5, fill=True, edgecolor='none',facecolor=colour_dict[var.upper()])
            else:
                rect = patches.Rectangle((left_of_box,bottom_of_box), spacing*0.8,  y_inc,alpha=0.5, fill=True, edgecolor='none',facecolor="dimgrey")

            ax.add_patch(rect)

            # sequence variant text
            ax.text(position, y_pos*y_inc, var, size=9, ha="center", va="center")

        # reference variant text

        ax.text(position, y_inc * -0.2, ref, size=9, ha="center", va="center")

        #polygon showing mapping from genome to spaced out snps
        x = [snp-0.5,snp+0.5,right_of_box,left_of_box,snp-0.5]

        y = [bottom_polygon,bottom_polygon,top_polygon,top_polygon,bottom_polygon]
        coords = list(zip(x, y))

        # draw polygon
        poly = patches.Polygon(coords, alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
        ax.add_patch(poly)

        rect = patches.Rectangle((left_of_box,top_polygon), spacing*0.8, y_inc,alpha=0.1, fill=True, edgecolor='none',facecolor="dimgrey")
        ax.add_patch(rect)

    if len(snp_dict) == 0:
        # snp position labels
        left_of_box = position-(0.4*position)
        right_of_box = position+(0.4*position)

        top_polygon = y_inc * -0.7
        bottom_polygon = y_inc * -1.7


    # reference variant rectangle
    rect = patches.Rectangle((0,(top_polygon)), length, y_inc ,alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
    ax.add_patch(rect)

    ax.text(-20,  y_inc * -0.2, label_map["reference"], size=9, ha="right", va="center")

    ref_genome_position = y_inc*-2.7

    # reference genome rectangle
    rect = patches.Rectangle((0,ref_genome_position), length, y_inc ,alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
    ax.add_patch(rect)

    for var in ref_vars:
        ax.plot([var,var],[ref_genome_position,ref_genome_position+(y_inc*0.98)], color="#cbaca4")


    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.yticks([])

    ax.set_xlim(0,length)
    if not flip_vertical:
        ax.set_ylim(ref_genome_position,y_level+(y_inc*1.05))
    else:
        ax.set_ylim(ref_genome_position,y_level+(y_inc*2.05))
        ax.invert_yaxis() # must be called after axis limits are set

    ax.tick_params(axis='x', labelsize=8)
    plt.xlabel("Position (base)", fontsize=9)
    plt.tight_layout()
    if not solid_background:
        plt.savefig(output, transparent=True)
    else:
        plt.savefig(output)

def get_colours(colour_palette):

    palettes = {"classic": {"A":"steelblue","C":"indianred","T":"darkseagreen","G":"skyblue"},
                "wes": {"A":"#CC8B3C","C":"#456355","T":"#541F12","G":"#B62A3D"},
                "primary": {"A":"green","C":"goldenrod","T":"steelblue","G":"indianred"},
                "purine-pyrimidine":{"A":"indianred","C":"teal","T":"teal","G":"indianred"},
                "greyscale":{"A":"#CCCCCC","C":"#999999","T":"#666666","G":"#333333"},
                "blues":{"A":"#3DB19D","C":"#76C5BF","T":"#423761","G":"steelblue"},
                "verity":{"A":"#EC799A","C":"#df6eb7","T":"#FF0080","G":"#9F0251"},
                "recombi":{"lineage_1":"steelblue","lineage_2":"#EA5463","Both":"darkseagreen","Private":"goldenrod"}
                }
    if colour_palette not in palettes:
        sys.stderr.write(red(f"Error: please select one of {palettes} for --colour-palette option\n"))
        sys.exit(-1)
    else:
        colour_dict = palettes[colour_palette]

    return colour_dict

def check_size_option(s):
    size_options = ["expand", "scale"]
    s_string = "\n - ".join(size_options)
    if s not in size_options:
        sys.stderr.write(red(f"Error: size option specified not one of:\n - {s_string}\n"))
        sys.exit(-1)

def check_format(f):
    formats = ["png", "jpg", "pdf", "svg", "tiff"]
    f_string = "\n - ".join(formats)
    if f not in formats:
        sys.stderr.write(red(f"Error: format specified not one of:\n - {f_string}\n"))
        sys.exit(-1)

def colour(text, text_colour):
    bold_text = 'bold' in text_colour
    text_colour = text_colour.replace('bold', '')
    underline_text = 'underline' in text_colour
    text_colour = text_colour.replace('underline', '')
    text_colour = text_colour.replace('_', '')
    text_colour = text_colour.replace(' ', '')
    text_colour = text_colour.lower()
    if 'red' in text_colour:
        coloured_text = RED
    elif 'green' in text_colour:
        coloured_text = GREEN
    elif 'yellow' in text_colour:
        coloured_text = YELLOW
    elif 'dim' in text_colour:
        coloured_text = DIM
    elif 'cyan' in text_colour:
        coloured_text = 'cyan'
    else:
        coloured_text = ''
    if bold_text:
        coloured_text += BOLD
    if underline_text:
        coloured_text += UNDERLINE
    if not coloured_text:
        return text
    coloured_text += text + END_FORMATTING
    return coloured_text

def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING
