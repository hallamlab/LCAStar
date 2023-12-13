#!python
"""
subsample_ncbi.py

Created by Niels Hanson on 2013-08-16.
Copyright (c) 2013 Steven J. Hallam Laboratory. All rights reserved.
"""

from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright 2013"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Niels W Hanson"
__status__ = "Release"


try:
     import os
     import re
     import argparse
     import glob
     import random
     from os import makedirs, sys, remove
     from sys import path
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)

what_i_do = "Create a contig->map file"
parser = argparse.ArgumentParser(description=what_i_do)

parser.add_argument('-i', dest='mp_mapping_file', type=str, nargs='?',
             required=True, help='file from preprocessed directory', default=None)                    
parser.add_argument('-n', dest='ncbi_summary', type=str, nargs='?',
             required=True, help='ncbi summary.txt file', default=None)
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
             required=True, help='output file of simulated metagenome', default=None)


def strip_decimal(text):
    return re.sub("\.[0-9]*","",text)

def main(argv):
    args = vars(parser.parse_args())
    
    # contig_to_ncbi_id
    contig_to_ncbi_id = {}
    
    # read in contig_map file
    with open(args['mp_mapping_file'],"r") as map_fh:
        for line in map_fh:
            fields = line.split("\t")
            fields = map(str.strip, fields)
            fields = map(str.strip, fields, "\n")
            if fields[0] not in contig_to_ncbi_id:
                ncbi_id = strip_decimal(fields[1])
                contig_to_ncbi_id[fields[0]] = ncbi_id
    
    ncbi_id_to_taxa = {}
    
    with open(args['ncbi_summary']) as ncbi_file:
        for line in ncbi_file:
            fields = line.split("\t")
            fields = map(str.strip, fields)
            fields = map(str.strip, fields, "\n")
            ncbi_id = strip_decimal(fields[0])
            if ncbi_id not in ncbi_id_to_taxa:
                ncbi_id_to_taxa[ncbi_id] = fields[5]
    
    # calculate and writeout resluts
    
    output = open(args['output_file'], 'w')
    
    output.write("Contig\tTaxonomy\n")
    
    for contig in contig_to_ncbi_id:
        fields = contig_to_ncbi_id[contig].split("|")
        if fields[3] in ncbi_id_to_taxa:
            line = "\t".join([contig, ncbi_id_to_taxa[fields[3]]]) + "\n"
            output.write(line)
    
    output.close()
    
    exit()  

# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
