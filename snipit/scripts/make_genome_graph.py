#!/usr/bin/env python3
import os
import argparse
import matplotlib as mpl
from matplotlib import pyplot as plt
import csv
import collections
import matplotlib.patches as patches
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from itertools import cycle
import math

colour_list = ["lightgrey","white"]
colour_dict = {"A":"steelblue","C":"indianred","T":"darkseagreen","G":"skyblue"}
colour_cycle = cycle(colour_list)

def parse_args():
    parser = argparse.ArgumentParser(description='Find sequences relative to Wuhan4 reference.')

    parser.add_argument("--input", action="store", type=str, dest="input")
    parser.add_argument("--ambiguities", action="store", type=str, dest="ambiguities")
    parser.add_argument("--output", action="store", type=str, dest="output")

    return parser.parse_args()

def next_colour():
    return next(colour_cycle)


def make_graph():

    args = parse_args()

    num_seqs = 0
    with open(args.input,newline="") as f:
        reader=csv.DictReader(f, delimiter='\t')
        for row in reader:
            num_seqs +=1

    height = math.sqrt(num_seqs)*2
    fig, ax = plt.subplots(1,1, figsize=(12,height), dpi=250)

    y_position = 0
    ref_vars = {}
    snp_dict = collections.defaultdict(list)

    
    amb_dict = collections.defaultdict(list)
    with open(args.ambiguities, newline="") as famb:
        reader = csv.DictReader(famb, delimiter="\t")
        for row in reader:
            if not row["ambiguous_snps"] == "":
                ambs = row["ambiguous_snps"].split(";")
                print(ambs)
                for ambiguity in ambs:
                    amb_dict[row["name"]].append(ambiguity)

    with open(args.input,newline="") as f:
        reader=csv.DictReader(f, delimiter='\t')
        for row in reader:
            y_position +=1
            snps = row["snps"].split(";")
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
                snp_dict[x_position].append((row["name"], ref, base, y_position))
                # ax.text(x_position, y_position, base, size=8, ha="center", va="center")
            if row["name"] in amb_dict:
                for amb in amb_dict[row["name"]]:
                    x_position = int(amb[:-2])
                    base = amb[-1]
                    ref = amb[-2]
                    ref_vars[x_position]=ref
                    snp_dict[x_position].append((row["name"], ref, base, y_position))

            ax.text(-20, y_position, row["name"], size=9, ha="right", va="center")
    
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
    plt.savefig(args.output)

if __name__ == '__main__':

    make_graph()