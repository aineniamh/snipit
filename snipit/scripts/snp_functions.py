#!/usr/bin/env python3

# imports of built-ins
import os
import sys
import argparse
import collections
from itertools import cycle
import csv
import math

# imports from other modules
from Bio import SeqIO
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

def qc_alignment(alignment,reference,cwd):
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
        if reference.split(".")[-1] not in ["gb","genbank"]:
            sys.stderr.write(red(f"Error: alignment file must contain more than just the reference. Either provide a reference genbank file or add more sequences to your alignment.\n"))
            sys.exit(-1)

    if len(set(lengths))!= 1:
        sys.stderr.write(red("Error: not all of the sequences in the alignment are the same length\n"))
        for i in lengths_info:
            print(f"{i[0]}\t{i[1]}\n")
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

def find_snps(reference_seq,input_seqs):

    non_amb = ["A","T","G","C"]
    snp_dict = {}
    record_snps = {}
    snp_counter = collections.Counter()
    for query_seq in input_seqs:
        snps =[]

        for i in range(len(query_seq)):
            bases = [query_seq[i],reference_seq[i]]
            if bases[0] != bases[1]:
                if bases[0] in non_amb and bases[1] in non_amb:
                    
                    snp = f"{i+1}{bases[1]}{bases[0]}" # position-reference-query
                    snp_counter[snp]+=1
                    snps.append(snp)
        snp_dict[query_seq] = snps

        for record in input_seqs[query_seq]:
            record_snps[record] = snps

    return snp_dict,record_snps,len(snp_counter)

def find_ambiguities(alignment, snp_dict):

    snp_sites = collections.defaultdict(list)
    for seq in snp_dict:
        snps = snp_dict[seq]
        for snp in snps:
            x_position = int(snp[:-2])-1
            ref = snp[-2]
            snp_sites[x_position]=ref

    amb_dict = {}
    
    for query_seq in alignment:
        snps =[]

        for i in snp_sites:
            bases = [query_seq[i],snp_sites[i]]
            if bases[0] != bases[1]:
                if bases[0] not in ["A","T","G","C"]:
                    
                    snp = f"{i+1}{bases[1]}{bases[0]}" # position-outgroup-query
                    snps.append(snp)
        
        for record in alignment[query_seq]:
            amb_dict[record] = snps

    return amb_dict

def make_graph(num_seqs,num_snps,amb_dict,snp_records,output,label_map,colour_dict,length,width,height,size_option):
    
    y_level = 0
    ref_vars = {}
    snp_dict = collections.defaultdict(list)
    
    for record in snp_records:

        # y level increments per record
        y_level +=1

        # for each record get the list of snps
        snps = snp_records[record]
        x = []
        y = []
        
        for snp in snps:
            # snp => 12345AT
            x_position = int(snp[:-2])
            base = snp[-1]
            ref = snp[-2]
            ref_vars[x_position]=ref
            # Add name of record, ref, SNP in record, y_level
            snp_dict[x_position].append((record, ref, base, y_level))

        # if there are ambiguities in that record, add them to the snp dict too
        if record in amb_dict:
            for amb in amb_dict[record]:
                # amb => 12345AN
                x_position = int(amb[:-2])
                base = amb[-1]
                ref = amb[-2]
                ref_vars[x_position]=ref
                # Add name of record, ref, SNP in record, y_level
                snp_dict[x_position].append((record, ref, base, y_level))

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
    # width and height of the figure
    fig, ax = plt.subplots(1,1, figsize=(width,height), dpi=250)

    y_level =0

    for record in snp_records:

        # y position increments
        y_level += y_inc

        # either grey or white
        col = next_colour()

        # for each record (sequence) draw a rectangle the length of the whole genome (either grey or white)
        rect = patches.Rectangle((0,y_level-(0.5*y_inc)), length, y_inc ,alpha=0.3, fill=True, edgecolor='none',facecolor=col)
        ax.add_patch(rect)

        # for each record add the name to the left hand side
        ax.text(-20, y_level, label_map[record], size=9, ha="right", va="center")

    position = 0
    for snp in sorted(snp_dict):
        position += spacing
        ax.text(position, y_level+(0.5*y_inc), snp, size=9, ha="center", va="bottom", rotation=90)
        # snp position labels
        left_of_box = position-(0.4*spacing)
        right_of_box = position+(0.4*spacing)

        top_polygon = y_inc * -0.7
        bottom_polygon = y_inc * -1.7

        for sequence in snp_dict[snp]:
            
            name,ref,var,y_pos = sequence
            bottom_of_box = (y_pos*y_inc)-(0.5*y_inc)
            

            # draw box for snp
            if var in colour_dict:
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

        # 
        rect = patches.Rectangle((left_of_box,top_polygon), spacing*0.8, y_inc,alpha=0.1, fill=True, edgecolor='none',facecolor="dimgrey")
        ax.add_patch(rect)


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
    ax.set_ylim(ref_genome_position,y_level+(y_inc*1.05))
    ax.tick_params(axis='x', labelsize=8)
    plt.xlabel("Genome position (base)", fontsize=9)
    plt.tight_layout()
    plt.savefig(output)

def get_colours(colour_palette):
    
    palettes = {"classic": {"A":"steelblue","C":"indianred","T":"darkseagreen","G":"skyblue"},
                "wes": {"A":"#CC8B3C","C":"#456355","T":"#541F12","G":"#B62A3D"}, 
                "primary": {"A":"green","C":"goldenrod","T":"steelblue","G":"indianred"},
                "purine-pyrimidine":{"A":"indianred","C":"teal","T":"teal","G":"indianred"},
                "greyscale":{"A":"#CCCCCC","C":"#999999","T":"#666666","G":"#333333"},
                "blues":{"A":"#3DB19D","C":"#76C5BF","T":"#423761","G":"steelblue"},
                "verity":{"A":"#EC799A","C":"#df6eb7","T":"#FF0080","G":"#9F0251"}
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
