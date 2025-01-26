# snipit
Summarise snps relative to a reference sequence

<img src="./docs/genome_graph.png" width="700">

### Install

```
pip install snipit
```

### Example Usage
Link to test data: [test.fasta](./docs/test.fasta) and [aa_test.fasta](./docs/aa_alignment.fasta)

- Basic usage for nucleotide alignments:
```
snipit test.fasta \
--output-file test
```
Default format output is `png`. Only specify output path/name (not extension).

- To change output format, use `--format`:
```
snipit test.fasta \
--output-file test \
--format pdf
```
Options: `png`, `jpg`, `pdf`, `svg`, `tiff`.

- To change color scheme, use `--colour-palette`:
```
snipit test.fasta \
--output-file test \
--colour-palette classic_extended
```

Other colours schemes:
```
classic, classic_extended, primary, purine-pyrimidine, greyscale, wes,verity, ugene
```
Use `ugene` for protein (aa) alignments.
Use `classic_extended` for colouring ambiguous bases.

- There are multiple options to control which SNPs or indels are included/excluded: 
```
snipit test.fasta \
--show-indels \
--include-positions '100-150' \
--exclude-positions '223 224 225'
```

- For control over ambiguous bases, use `--ambig-mode` to specify how ambiguous bases are handled:
```
[all] include all ambig such as N,Y,B in all positions
[snps] only include ambig if a snp is present at the same position - Default 
[exclude] remove all ambig, same as depreciated --exclude-ambig-pos
```
Use the colour palette `classic_extended` when plotting with `all` or `snps`.

 - Recombination mode is designed to assist with recombination analysis for SC2. This mode allows for colouring of mutations present in two references. For recombination mode, three flags are required: `--reference`,`--recombi-mode`,`--recombi-references`.

The specified `--reference` must be different from the `--recombi-references`.
```
snipit test.fasta \
--reference USA_3 \
--recombi-mode \
--recombi-references "USA_1,USA_2"
```

For amino acid alignments, specify the sequence type as `aa`, use the colour palette `ugene`:
```
snipit test.prot.fasta \
--sequence-type aa \
--colour-palette ugene \
--output-file test.prot
```

There are several more options, see below for full usage.

### Issues

If you see an error like: 
```
ModuleNotFoundError: No module named 'pkg_resources'
```
This may mean your python install did not come with setuptools. Install [setuptools](https://packaging.python.org/en/latest/guides/installing-using-linux-tools/) and this should resolve the issue.

### Full Usage
```
snipit

optional arguments:
  -h, --help            show this help message and exit

Input options:
  alignment             Input alignment fasta file
  -t {nt,aa}, --sequence-type {nt,aa}
                        Input sequence type: aa or nt
  -r REFERENCE, --reference REFERENCE
                        Indicates which sequence in the alignment is the
                        reference (by sequence ID). Default: first sequence in
                        alignment
  -l LABELS, --labels LABELS
                        Optional csv file of labels to show in output snipit
                        plot. Default: sequence names
  --l-header LABEL_HEADERS
                        Comma separated string of column headers in label csv.
                        First field indicates sequence name column, second the
                        label column. Default: 'name,label'

Mode options:
  --recombi-mode        Allow colouring of query seqeunces by mutations
                        present in two 'recombi-references' from the input
                        alignment fasta file
  --recombi-references RECOMBI_REFERENCES
                        Specify two comma separated sequence IDs in the input
                        alignment to use as 'recombi-references'. Ex.
                        Sequence_ID_A,Sequence_ID_B
  --cds-mode            Assumes sequence supplied is a coding sequence

Output options:
  -d OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory. Default: current working directory
  -o OUTFILE, --output-file OUTFILE
                        Output file name stem. Default: snp_plot
  -s, --write-snps      Write out the SNPs in a csv file.
  -f FORMAT, --format FORMAT
                        Format options (png, jpg, pdf, svg, tiff) Default: png

Figure options:
  --height HEIGHT       Overwrite the default figure height
  --width WIDTH         Overwrite the default figure width
  --size-option SIZE_OPTION
                        Specify options for sizing. Options: expand, scale
  --solid-background    Force the plot to have a solid background, rather than
                        a transparent one.
  -c , --colour-palette 
                        Specify colour palette. Options: [classic,
                        classic_extended, primary, purine-pyrimidine,
                        greyscale, wes, verity, ugene]. Use ugene for protein
                        alignments.
  --flip-vertical       Flip the orientation of the plot so sequences are
                        below the reference rather than above it.
  --sort-by-mutation-number
                        Render the graph with sequences sorted by the number
                        of SNPs relative to the reference (fewest to most).
                        Default: False
  --sort-by-id          Sort sequences alphabetically by sequence id. Default:
                        False
  --sort-by-mutations SORT_BY_MUTATIONS
                        Sort sequences by bases at specified positions.
                        Positions are comma separated integers. Ex. '1,2,3'
  --high-to-low         If sorted by mutation number is selected, show the
                        sequences with the fewest SNPs closest to the
                        reference. Default: False

SNP options:
  --show-indels         Include insertion and deletion mutations in snipit
                        plot.
  --include-positions INCLUDED_POSITIONS [INCLUDED_POSITIONS ...]
                        One or more range (closed, inclusive; one-indexed) or
                        specific position only included in the output. Ex.
                        '100-150' or Ex. '100 101' Considered before '--
                        exclude-positions'.
  --exclude-positions EXCLUDED_POSITIONS [EXCLUDED_POSITIONS ...]
                        One or more range (closed, inclusive; one-indexed) or
                        specific position to exclude in the output. Ex.
                        '100-150' or Ex. '100 101' Considered after '--
                        include-positions'.
  --ambig-mode {all,snps,exclude}
                        Controls how ambiguous bases are handled - [all]
                        include all ambig such as N,Y,B in all positions;
                        [snps] only include ambig if a snp is present at the
                        same position; [exclude] remove all ambig, same as
                        depreciated --exclude-ambig-pos

Misc options:
  -v, --version         show program's version number and exit
```

### Cite

Please cite this tool as follows:
```
Aine O'Toole, snipit (2024) GitHub repository, https://github.com/aineniamh/snipit
```
