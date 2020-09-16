# snipit
Summarise snps relative to a reference sequence


<img src="./docs/genome_graph.png" width="700">

### Usage
```
usage: snipit <alignment> [options]

snipit 


positional arguments:
  alignment                 Input alignment fasta file

optional arguments:
  -h, --help                show this help message and exit
  -r REFERENCE, 
  --reference REFERENCE
                            Indicates which sequence in the alignment is the reference (by sequence ID). Can also accept a genbank file. Default: first sequence in alignment
  -o OUTPUT_DIR, 
  --output-dir OUTPUT_DIR
                            Output directory. Default: current working directory
  -f FORMAT, 
  --format FORMAT
                            Format options. Default: png
  --height HEIGHT           Overwrite the default figure height
  --width WIDTH             Overwrite the default figure width
```