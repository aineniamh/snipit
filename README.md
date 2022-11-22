# snipit
Summarise snps relative to a reference sequence


<img src="./docs/genome_graph.png" width="700">

### Usage
```
usage: snipit <alignment> [options]

snipit

positional arguments:
  alignment             Input alignment fasta file

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Indicates which sequence in the alignment is the reference (by sequence ID). Default: first sequence in
                        alignment
  -l LABELS, --labels LABELS
                        Optional csv file of labels to show in output snipit plot. Default: sequence names
  --l-header LABEL_HEADERS
                        Comma separated string of column headers in label csv. First field indicates sequence name column, second
                        the label column. Default: 'name,label'
  -d OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory. Default: current working directory
  -o OUTFILE, --output-file OUTFILE
                        Output file name stem. Default: snp_plot
  -s, --write-snps      Write out the SNPs in a csv file.
  -f FORMAT, --format FORMAT
                        Format options (png, jpg, pdf, svg, tiff) Default: png
  --height HEIGHT       Overwrite the default figure height
  --width WIDTH         Overwrite the default figure width
  --size-option SIZE_OPTION
                        Specify options for sizing. Options: expand, scale
  --solid-background    Force the plot to have a solid background, rather than a 
                        transparent one.
  --flip-vertical       Flip the orientation of the plot so sequences are below the 
                        reference rather than above it.
  --snps-only           Ignore insertion and deletion mutations and only plot SNPs
                        (legacy behaviour).
  --include-positions INCLUDED_POSITIONS [INCLUDED_POSITIONS ...]
                        One or more range (closed, inclusive; one-indexed) or specific position only included in the output. Ex.
                        '100-150' or Ex. '100 101' Considered before '--exclude-positions'.
  --exclude-positions EXCLUDED_POSITIONS [EXCLUDED_POSITIONS ...]
                        One or more range (closed, inclusive; one-indexed) or specific position to exclude in the output. Ex.
                        '100-150' or Ex. '100 101' Considered after '--include-positions'.
  --exclude-ambig-pos   Exclude positions with ambig base in any sequences. Considered 
                        after '--include-positions'
  --sort-by-mutation-number
                        Render the graph with sequences sorted by the number of SNPs relative to the reference (fewest to most).
                        Default: False
  --sort-by-id          Sort sequences alphabetically by sequence id. Default: False
  --sort-by-mutations SORT_BY_MUTATIONS
                        Sort sequences by bases at specified positions. Positions are comma separated integers. Ex. '1,2,3'
  --high-to-low         If sorted by mutation number is selected, show the sequences 
                        with the fewest SNPs closest to the
                        reference. Default: False
  -v, --version         show program's version number and exit
  -c COLOUR_PALETTE, --colour-palette COLOUR_PALETTE
                        Specify colour palette. Options: primary, classic, purine-pyrimidine, greyscale, wes, verity
  --recombi-mode        Allow colouring of query seqeunces by mutations present in two 
                        'recombi-references' from the input
                        alignment fasta file
  --recombi-references RECOMBI_REFERENCES
                        Specify two comma separated sequence IDs in the input alignment to use as 'recombi-references'. Ex.
                        Sequence_ID_A,Sequence_ID_B
```

### Install

```
pip install snipit
```
