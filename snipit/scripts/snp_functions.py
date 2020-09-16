#!/usr/bin/env python3
import os
import argparse
import matplotlib as mpl
from matplotlib import pyplot as plt
import csv
import collections
import matplotlib.patches as patches
from matplotlib.patches import Polygon
from itertools import cycle
import math
import sys
from Bio import SeqIO

colour_list = ["lightgrey","white"]
colour_dict = {"A":"steelblue","C":"indianred","T":"darkseagreen","G":"skyblue"}
colour_cycle = cycle(colour_list)



def qc_alignment(alignment):
    lengths = []
    lengths_info = []
    num_seqs = 0
    
    record_ids = []
    ref_input = ""
    for record in SeqIO.parse(alignment, "fasta"):
        if ref_input == "":
            ref_input = record.id
        lengths.append(len(record))
        record_ids.append(record.id)
        lengths_info.append((record.id, len(record)))
        num_seqs +=1

    if len(set(lengths))!= 1:
        sys.stderr.write("Error: not all of the sequences in the alignment are the same length\n")
        for i in lengths_info:
            print(f"{i[0]}\t{i[1]}\n")
        sys.exit(-1)
    
    return num_seqs,ref_input,record_ids

def reference_qc(reference, record_ids,cwd):
    ref_file = ""
    if "." in reference and reference.split(".")[-1] in ["gb","genbank"]:
        ref_path = os.path.join(cwd, reference)
        if not os.path.exists(ref_path):
            sys.stderr.write(f"Error: can't find genbank file at {ref_path}\n")
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
                sys.stderr.write(f"Error: more than one record found in reference genbank file\n")
                sys.exit(-1)

    elif reference not in record_ids:
        sys.stderr.write(f"Error: input reference {reference} not found in alignment\n")
        sys.exit(-1)

    else:
        ref_input = reference

    return ref_file, ref_input

def next_colour():
    return next(colour_cycle)

def get_ref_and_alignment(input_file,reference):
    input_seqs = collections.defaultdict(list)
    reference_seq = ""

    for record in SeqIO.parse(input_file, "fasta"):
        if record.id == reference:
            reference_seq = record.seq.upper()
        else:
            input_seqs[str(record.seq).upper()].append(record.id)

    return reference_seq, input_seqs

def find_snps(reference_seq,input_seqs):

    non_amb = ["A","T","G","C"]
    snp_dict = {}
    record_snps = {}

    for query_seq in input_seqs:
        snps =[]

        for i in range(len(query_seq)):
            bases = [query_seq[i],reference_seq[i]]
            if bases[0] != bases[1]:
                if bases[0] in non_amb and bases[1] in non_amb:
                    
                    snp = f"{i+1}{bases[1]}{bases[0]}" # position-reference-query
                    snps.append(snp)
        snp_dict[query_seq] = snps

        for record in input_seqs[query_seq]:
            record_snps[record] = snps

    return snp_dict,record_snps

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

def make_graph(num_seqs,input_file,amb_dict,snp_records,output):

    height = math.sqrt(num_seqs)*2
    fig, ax = plt.subplots(1,1, figsize=(12,height), dpi=250)

    y_position = 0
    ref_vars = {}
    snp_dict = collections.defaultdict(list)

    for record in snp_records:
        y_position +=1
        snps = snp_records[record]
        x = []
        y = []
        col = next_colour()
        rect = patches.Rectangle((0,y_position-0.5), 29903, 1 ,alpha=0.3, fill=True, edgecolor='none',facecolor=col)
        ax.add_patch(rect)
        
        for snp in snps:
            x_position = int(snp[:-2])
            base = snp[-1]
            ref = snp[-2]
            ref_vars[x_position]=ref
            snp_dict[x_position].append((record, ref, base, y_position))
            # ax.text(x_position, y_position, base, size=8, ha="center", va="center")
        if record in amb_dict:
            for amb in amb_dict[record]:
                x_position = int(amb[:-2])
                base = amb[-1]
                ref = amb[-2]
                ref_vars[x_position]=ref
                snp_dict[x_position].append((record, ref, base, y_position))

        ax.text(-20, y_position, record, size=9, ha="right", va="center")
    
    spacing = 29903/(len(snp_dict)+1)
    
    position = 0
    for snp in sorted(snp_dict):
        position += spacing
        ax.text(position, y_position+1, snp, size=9, ha="center", va="bottom", rotation=90)

        # snp position labels
        
        for sequence in snp_dict[snp]:
            # sequence variant text
            name,ref,var,y_pos = sequence
            if var in colour_dict:
                rect = patches.Rectangle((position-(0.4*spacing),y_pos-0.5), spacing*0.8, 1 ,alpha=0.5, fill=True, edgecolor='none',facecolor=colour_dict[var.upper()])
            else:
                rect = patches.Rectangle((position-(0.4*spacing),y_pos-0.5), spacing*0.8, 1 ,alpha=0.5, fill=True, edgecolor='none',facecolor="dimgrey")
            ax.add_patch(rect)
            ax.text(position, y_pos, var, size=9, ha="center", va="center")

        # reference variant text
        ax.text(position, -0.2, ref, size=9, ha="center", va="center") 

        #polygon showing mapping from genome to spaced out snps
        x = [snp-0.5,snp+0.5,position+(0.4*spacing),position-(0.4*spacing),snp-0.5]
        y = [-1.7,-1.7,-0.7,-0.7,-1.7]
        coords = list(zip(x, y))
        poly = patches.Polygon(coords, alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
        ax.add_patch(poly)
        rect = patches.Rectangle((position-(0.4*spacing),-0.7), spacing*0.8, 1 ,alpha=0.1, fill=True, edgecolor='none',facecolor="dimgrey")
        ax.add_patch(rect)

    # reference variant rectangle
    rect = patches.Rectangle((0,-0.7), 29903, 1 ,alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
    ax.add_patch(rect)
    ax.text(-20, -0.2, "Reference", size=9, ha="right", va="center")
    # reference genome rectangle
    rect = patches.Rectangle((0,-2.7), 29903, 1 ,alpha=0.2, fill=True, edgecolor='none',facecolor="dimgrey")
    ax.add_patch(rect)

    for var in ref_vars:
        ax.plot([var,var],[-2.69,-1.71], color="#cbaca4")


    ax.spines['top'].set_visible(False) ## make axes invisible
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)

    plt.yticks([])
            
    ax.set_xlim(0,29903)
    ax.set_ylim(-2.7,y_position+1)
    ax.tick_params(axis='x', labelsize=8)
    plt.xlabel("Genome position (base)", fontsize=9)
    plt.tight_layout()
    plt.savefig(output)

