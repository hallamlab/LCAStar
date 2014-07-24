#!python
"""
create_analysis_table.py

Created by Niels Hanson
Copyright (c) 2014 Steven J. Hallam Laboratory. All rights reserved.
"""

from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright 2014"
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
     from LCAStar import *
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)

what_i_do = "Creates LCA_Star analysis table from MetaPathways output for downstream analysis."
parser = argparse.ArgumentParser(description=what_i_do)

parser.add_argument('-i', dest='input_mapping', type=str, nargs='?',
             required=True, help='MetaPathways: fasta input mapping file from preprocessed/', default=None)                    
parser.add_argument('--ft_table', dest='ft_table', type=str, nargs='?',
             required=False, help='MetaPathways: functional and taxonoic table', default=None)
parser.add_argument('--sum_table', dest='sum_table', type=str, nargs='?',
             required=False, help='NCBI Genomes table: summary.txt', default=None)
parser.add_argument('--ncbi_tree', dest='ncbi_tree', type=str, nargs='?',
             required=False, help='MetaPathways: NCBI tree file', default=None)
parser.add_argument('-o', dest='output', type=str, nargs='?',
             required=False, help='output file', default=None)

# Functional classes    
"""
A class to represent an individual .fasta record within a fasta file.
It contains three fields: (1) longname - representing the whole annotation 
contained after the carat symbol (>), (2) sequence - the amino acid or nucleotide 
sequence, and (3) name a shortened name based on the information at appears before
first space in the annotation line.
"""
class FastaRecord():
 def __init__(self, longname, sequence):
   self.longname = longname
   self.sequence = sequence
   fields = [ x.strip() for x in self.longname.split(' ') ]
   if len(fields) > 0:
      self.name = fields[0]
   else:
      self.name = None

"""
Takes a fasta file as input and returns a list of FastaRecords for
each correctly parsed record in the file.
"""
class FastaReader():
 stop = False
 START_PATTERN = re.compile(r'^>')
 name = None
 future_name =None
 sequence=""
 def __init__(self, fasta_filename):
     try:
         self.file = open(fasta_filename, 'r')
         # print fasta_filename
     except IOError:
         print "Cannot open fasta file " + fasta_filename

 def __iter__(self):
     return self


 def next(self):
     if self.stop:
       raise StopIteration

     try:
        if not self.name: 
            self.name = self.file.readline().strip()
        line = self.file.readline().strip()
     except:
        line = None


     if not line:
        self.stop = True
        raise StopIteration


     fragments = []
     while line and not self.START_PATTERN.search(line):
         fragments.append(line.strip()) 
         line = self.file.readline()

    # print line
     if self.future_name:
         self.name = self.future_name

     if line:
       self.future_name = line.strip()

     self.sequence =''.join(fragments)
     self.seqname = self.name

     return FastaRecord(self.name, self.sequence)

def main(argv):
    args = vars(parser.parse_args())
    
    # read the input fasta_files to be mapped to create read2origin map
    fh = open(args["input_mapping"], "r")
    lines = fh.readlines()
    fh.close()

    read2origin = {} # hash mapping input sequences to their original header
    
    # fill read2origin
    for l in lines:
       fields =l.split("\t")
       fields = map(str.strip, fields)
       read_name = fields[0]
       gi_line = fields[1].strip("\n")
       gi_line = gi_line.split("|")
       ncbi_id = re.sub("\.[0-9]+$", "", gi_line[3] )
       read2origin[fields[0]] = ncbi_id
    
    # Parse NCBI: summary table
    fh = open(args["sum_table"], "r")
    lines = fh.readlines()
    fh.close()
    
    ncbiID_to_taxaID = {} # ncbiID to NCBI Tree taxaID
    taxaID_to_taxa = {} # ncbiID to full taxa name
    for l in lines:
       fields = l.split("\t")
       fields = map(str.strip, fields)
       ncbi_id = re.sub("\.[0-9]+$", "", fields[0])
       tax_id = fields[3]
       tax = fields[5]
       if ncbi_id not in ncbiID_to_taxaID:
          ncbiID_to_taxaID[ncbi_id] = tax_id
       if tax_id not in taxaID_to_taxa:
          taxaID_to_taxa[tax_id] = tax
    
    # read functional and taxonomic table
    fh = open(args["ft_table"], "r")
    header = fh.readline()
    header = header.split("\t") # list of headers
    lines = fh.readlines()

    contig_to_orfs = {} # get a list of ORFs for a specific contig
    orfs_to_lca = {} # get the lca taxonomy

    for l in lines:
       fields = l.split("\t")
       # ORF_ID   ORF_length  start   end Contig_Name Contig_length   strand  ec  taxonomy    product
       orf_id = fields[0]
       contig = fields[4]
       lca = fields[8]
       if contig not in contig_to_orfs:
          contig_to_orfs[contig] = []
       contig_to_orfs[contig].append(orf_id)
       if orf_id not in orfs_to_lca:
          orfs_to_lca[orf_id] = lca
       # get original taxonomy
       # print taxaID_to_taxa[ncbiID_to_taxaID[read2origin[contig]]]
    
    # Build the LCA Star NCBI Tree
    print "Loading LCAStar:"
    lcastar = LCAStar(args["ncbi_tree"])
    print "Done."
    
    # set LCAStar parameters
    lcastar.setLCAStarParameters(min_depth = 2, alpha = 0.5, min_reads = 3 )
    
    # small helper to translate the list
    def get_orfs_taxa(orfs):
        list = []
        for o in orfs:
            list.append(orfs_to_lca[o])
        return list
    output = open(args["output"], "w") 
    header = "\t".join(["contig","real","taxa","method","dist", "wtd", "real_linage"])
    output.write(header + "\n")
    for c in contig_to_orfs:
        # the read taxa
        contig = c
        real = taxaID_to_taxa[ncbiID_to_taxaID[read2origin[c]]]
        taxa_list = get_orfs_taxa(contig_to_orfs[c])
        lca_list = []
        for t in taxa_list:
             lca_list.append([ t ])
        taxon = lcastar.lca_star(taxa_list)
        real_lineage = lcastar.get_lineage(lcastar.get_a_Valid_ID([real]))
        real_lineage = map( str, map(lcastar.translateIdToName, real_lineage))
        line = "\t".join([contig, real, taxon, "LCA_Star", str(lcastar.get_distance(taxon, real)), str(lcastar.wtd_distance(real, taxon))] + real_lineage[::-1])
        output.write(line + "\n")
        taxon = lcastar.lca_majority(taxa_list) 
        line = "\t".join([contig, real, taxon, "Majority", str(lcastar.get_distance(taxon, real)), str(lcastar.wtd_distance(real, taxon))] + real_lineage[::-1])
        output.write(line + "\n")
        taxon = lcastar.getTaxonomy(lca_list) 
        line = "\t".join([contig, real, taxon, "LCA_Squared", str(lcastar.get_distance(taxon, real)), str(lcastar.wtd_distance(real, taxon))] + real_lineage[::-1])
        output.write(line + "\n")

    output.close()
# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])