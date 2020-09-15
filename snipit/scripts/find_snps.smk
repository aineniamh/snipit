from Bio import Phylo
from Bio import SeqIO
import csv
import collections


config["tree_stems"] = config["catchment_str"].split(",")

rule all:
    input:
        expand(os.path.join(config["outdir"],"snp_reports","{tree}.snps.txt"),tree=config["tree_stems"]),
        expand(os.path.join(config["outdir"],"figures","genome_graph_{tree}.png"), tree=config["tree_stems"]),
        os.path.join(config["outdir"],"gather_prompt.txt")

rule extract_taxa:
    input:
        collapsed_tree = os.path.join(config["outdir"],"local_trees","{tree}.tree")
    output:
        tree_taxa = os.path.join(config["tempdir"], "local_trees_seqs","{tree}_taxon_names.txt")
    shell:
        "clusterfunk get_taxa -i {input.collapsed_tree:q} --in-format nexus -o {output.tree_taxa:q} --out-format newick"

rule gather_fasta_seqs:
    input:
        aligned_query_seqs = config["aligned_query_seqs"],
        cog_seqs = config["all_cog_seqs"],
        outgroup_fasta = config["outgroup_fasta"],
        tree_taxa = rules.extract_taxa.output.tree_taxa
    output:
        aln = os.path.join(config["tempdir"], "seqs_for_snps","{tree}.fasta")
    run:
        taxa = []
        with open(input.tree_taxa, "r") as f:
            for l in f:
                l = l.rstrip("\n")
                taxa.append(l)

        with open(output.aln, "w") as fw:
            for record in SeqIO.parse(input.outgroup_fasta, "fasta"):
                fw.write(f">outgroup\n{record.seq}\n")

            for record in SeqIO.parse(input.aligned_query_seqs, "fasta"):
                if record.id in taxa:
                    fw.write(f">{record.description}\n{record.seq}\n")

            for record in SeqIO.parse(input.cog_seqs,"fasta"):
                if record.id in taxa:
                    fw.write(f">{record.description}\n{record.seq}\n")


rule assess_snps:
    input:
        aln = rules.gather_fasta_seqs.output.aln
    params:
        tree = "{tree}"
    output:
        snp_report = os.path.join(config["outdir"], "snp_reports","{tree}.snps.txt")
    shell:
        """
        find_snps.py --input {input.aln:q} --output {output.snp_report:q} --tree {params.tree}
        """

rule ambiguities_at_snp_sites:
    input:
        seqs = os.path.join(config["tempdir"], "seqs_for_snps", "{tree}.fasta"),
        report = os.path.join(config["outdir"],"snp_reports", "{tree}.snps.txt")
    output:
        snp_report = os.path.join(config["outdir"], "snp_reports", "ambiguities_{tree}.snps.txt")
    shell:
        """
        find_ambiguities.py --input {input.seqs:q} --output {output.snp_report:q} --report {input.report:q}
        """

rule make_snp_figure:
    input:
        ambs = os.path.join(config["outdir"], "snp_reports", "ambiguities_{tree}.snps.txt"),
        snps = os.path.join(config["outdir"],"snp_reports","{tree}.snps.txt")
    output:
        os.path.join(config["outdir"],"figures","genome_graph_{tree}.png")
    shell:
        """
        make_genome_graph.py --input {input.snps:q} --ambiguities {input.ambs:q} --output {output[0]} 
        """

rule gather_graphs:
    input:
        expand(os.path.join(config["outdir"],"figures","genome_graph_{tree}.png"), tree=config["tree_stems"])
    output:
        os.path.join(config["outdir"],"gather_prompt.txt")
    shell:
        "touch {output}"
    