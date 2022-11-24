#!/usr/bin/env python3

# imports of built-ins

import os
import sys

from snipit.utils.colours import *

def check_ref(recombi_mode):
    if recombi_mode:
        sys.stderr.write(red(f"Error: Please explicitly state reference sequence when using `--recombi-mode`\n"))
        sys.exit(-1)
        

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

def amend_order(record_order,recombi_refs):
    
    # Reorder list to put recombi_references at the start
    record_order.remove(recombi_refs[0])
    record_order.insert(0, recombi_refs[0])
    record_order.remove(recombi_refs[1])
    record_order.insert(1, recombi_refs[1])

    return record_order

def run_recombi(record_order,recombi_references,snp_records):
    # Get a list of SNPs present in each recombi_reference
    recombi_snps,recombi_refs = recombi.recombi_ref_snps(recombi_references, snp_records)
    # Set the colour palette to "recombi"
    colour_dict = get_colours("recombi")

    record_order = amend_order(record_order,recombi_refs)