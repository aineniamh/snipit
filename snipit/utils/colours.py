#!/usr/bin/env python3

ND_FORMATTING = '\033[0m'
BOLD = '\033[1m'
UNDERLINE = '\033[4m'
RED = '\033[31m'
GREEN = '\033[32m'
YELLOW = '\033[93m'
CYAN = '\u001b[36m'
DIM = '\033[2m'

from itertools import cycle, chain


COLOUR_LIST = ["lightgrey","white"]
COLOUR_CYCLE = cycle(COLOUR_LIST)


def red(text):
    return RED + text + END_FORMATTING

def cyan(text):
    return CYAN + text + END_FORMATTING

def green(text):
    return GREEN + text + END_FORMATTING

def yellow(text):
    return YELLOW + text + END_FORMATTING


def next_colour():
    return next(COLOUR_CYCLE)


def get_colours(colour_palette):
    
    palettes = {"classic": {"A":"steelblue","C":"indianred","T":"darkseagreen","G":"skyblue"},
                "wes": {"A":"#CC8B3C","C":"#456355","T":"#541F12","G":"#B62A3D"}, 
                "primary": {"A":"green","C":"goldenrod","T":"steelblue","G":"indianred"},
                "purine-pyrimidine":{"A":"indianred","C":"teal","T":"teal","G":"indianred"},
                "greyscale":{"A":"#CCCCCC","C":"#999999","T":"#666666","G":"#333333"},
                "blues":{"A":"#3DB19D","C":"#76C5BF","T":"#423761","G":"steelblue"},
                "verity":{"A":"#EC799A","C":"#df6eb7","T":"#FF0080","G":"#9F0251"},
                "recombi":{"lineage_1":"steelblue","lineage_2":"#EA5463","Both":"darkseagreen","Private":"goldenrod"}
                }
    if colour_palette not in palettes:
        sys.stderr.write(red(f"Error: please select one of {palettes} for --colour-palette option\n"))
        sys.exit(-1)
    else:
        colour_dict = palettes[colour_palette]

    return colour_dict