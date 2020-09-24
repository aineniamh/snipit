# snipit
Summarise snps relative to a reference sequence


<img src="./docs/genome_graph.png" width="700">

### Usage
```
snipit

positional arguments:
  alignment             Input alignment fasta file

optional arguments:
  -h, --help            show this help message and exit
  -r REFERENCE, --reference REFERENCE
                        Indicates which sequence in the alignment is the reference (by sequence ID). Default: first sequence in alignment
  -l LABELS, --labels LABELS
                        Optional csv file of labels to show in output snipit plot. Default: sequence names
  --l-header LABEL_HEADERS
                        Comma separated string of column headers in label csv. First field indicates sequence name column, second the label column.
                        Default: 'name,label'
  -d OUTPUT_DIR, --output-dir OUTPUT_DIR
                        Output directory. Default: current working directory
  -o OUTFILE, --output-file OUTFILE
                        Output file name stem. Default: snp_plot
  -f FORMAT, --format FORMAT
                        Format options (png, jpg, pdf, svg, tiff) Default: png
  --height HEIGHT       Overwrite the default figure height
  --width WIDTH         Overwrite the default figure width
  --size-option SIZE_OPTION
                        Specify options for sizing. Options: expand, scale
  -c COLOUR_PALETTE, --colour-palette COLOUR_PALETTE
                        Specify colour palette. Options: primary, classic, purine-pyrimidine, greyscale, wes, verity
```

### Install

```
pip install git+https://github.com/aineniamh/snipit.git
```
