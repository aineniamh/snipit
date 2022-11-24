#!/usr/bin/env python3

# imports of built-ins

import os
import sys
from Bio import SeqIO

from snipit.utils.colours import *


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
        if reference:
            if reference.split(".")[-1] not in ["gb","genbank"]:
                sys.stderr.write(red(f"Error: alignment file must contain more than just the reference. Either provide a reference genbank file or add more sequences to your alignment.\n"))
                sys.exit(-1)
        else:
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
