#!python
"""
create_analysis_table.py

Created by Niels Hanson
Copyright (c) 2014 Steven J. Hallam Laboratory. All rights reserved.
"""

from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright (c) 2014 Steven J. Hallam Laboratory. All rights reserved."
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
     import functools
     from os import makedirs, sys, remove
     from sys import path
     from python_resources.LCAStar import *
     from python_resources.fastareader import *
except:
     print """ Could not load some modules """
     print """ """
     sys.exit(3)

what_i_do = """Computes three estimates of Contig Taxonomy: LCA^2, Majority, and entropy-based LCA* \
as described in the LCAStar paper (Hanson, et al. 2014).
"""
parser = argparse.ArgumentParser(description=what_i_do)

parser.add_argument('-i', '--parsed_blast_file', dest='parsed_blast', type=str, nargs='?', required=True,
    help='MetaPathways Output: Parsed (B)LAST annotation file.', default=None)
parser.add_argument('-m', '--mapping_file', dest='mapping_file', type=str, nargs='?', required=True,
    help='MetaPathways Output: Input mapping file in preprocessed directory.', default=None)
parser.add_argument('--ncbi_tree', dest='ncbi_tree', type=str, nargs='?', required=False,
    help='MetaPathways: NCBI tree file', default=None)
parser.add_argument('--ncbi_megan_map', dest='ncbi_megan_map', type=str, nargs='?', required=False,
    help='Preferred mapping for NCBI names (map.ncbi)', default=None)
parser.add_argument('-o', '--output', dest='output', type=str, nargs='?', required=False,
    help='output file of predicted taxonomies', default=None)

def translate_to_prefered_name(id, ncbi_megan_map, lcastar):
    id_str = str(id)
    if id_str in ncbi_megan_map:
        return ncbi_megan_map[id_str] + " (" + id_str + ")"
    else:
        res = lcastar.translateIdToName(id)
        if res:
            return res + " (" + id_str + ")"
        else:
            return "Unknown (" + id_str + ")"

# Takes a line and
def clean_tab_lines(line):
    fields = line.split("\t")
    fields = map(str.strip, fields)
    fields = map(str.strip, fields, "\n")
    return fields

def main(argv):
    args = vars(parser.parse_args())

    ## read the input fasta_files to be mapped to create read2origin map
    contig2origin = {} # hash mapping input sequences to their original header
    with open(args["mapping_file"], "r") as fh:
       for l in fh:
           fields = clean_tab_lines(l)
           read_name = fields[0]
           original_name = fields[1]
           contig2origin[read_name] = original_name

    ## create preferred mapping
    ncbi_megan_map = {} # hash map from given taxonomy to preferred one used by megan
    with open(args["ncbi_megan_map"], 'r') as meganfile:
        for line in meganfile:
             fields = line.split("\t")
             fields = map(str.strip, fields)
             ncbi_megan_map[fields[0]] = fields[1]


    # contig_taxa_datastructure
    contig_to_taxa = {}
    contig_pattern = re.compile("^(.*)_([0-9]*)$")
    taxonomy_pattern = re.compile("\[(.*)\]")

    # read blast table
    with open(args["parsed_blast"], "r") as fh:
        for l in fh:
            if re.match("^\#", l):
                clean_tab_lines(l)
            else:
                fields = clean_tab_lines(l)
                contig_hits = contig_pattern.search(fields[0])
                if contig_hits:
                    contig = contig_hits.group(1)
                    orf = contig_hits.group(2)
                    # add to data structre if it doesn't exist
                    if contig not in contig_to_taxa:
                        contig_to_taxa[contig] = {}
                    if orf not in contig_to_taxa[contig]:
                        contig_to_taxa[contig][orf] = []


                    # pull taxonomy out of annotation
                    taxa_hits = taxonomy_pattern.search(fields[9])
                    if taxa_hits:
                        taxa = taxa_hits.group(1)
                        contig_to_taxa[contig][orf].append([taxa])
                    else:
                        continue
                else:
                    continue



    ## Build the LCA Star NCBI Tree
    print "Loading LCAStar:"
    lcastar = LCAStar(args["ncbi_tree"])
    lcastar.setLCAStarParameters(min_depth = 1, alpha = 0.51, min_reads = 1)
    print "Done."

    ## Calculate LCA for each ORF
    contig_to_lca = {}
    for contig in contig_to_taxa:
        for orf in contig_to_taxa[contig]:
            if contig not in contig_to_lca:
                contig_to_lca[contig] = {}
            if orf not in contig_to_lca[contig]:
                contig_to_lca[contig][orf] = None
            #TODO: Think about a way to best-hit here as an alternative
            contig_to_lca[contig][orf] = lcastar.getTaxonomy(contig_to_taxa[contig][orf])

    ## calculate taxonomy statistics LCA,  for each ORF
    # contig_to_taxa = {}
    for contig in contig_to_lca:
        orf_lcas = []
        simple_list = []
        print contig
        for orf in contig_to_lca[contig]:
            lca = contig_to_lca[contig][orf]
            orf_lcas.append( [ lca ] )
            simple_list.append( lca )
            print "\t", orf, lca
        lca_squared = lcastar.getTaxonomy( orf_lcas )
        majority = lcastar.lca_majority( simple_list )
        lca_star = lcastar.lca_star( simple_list )
        print simple_list
        print contig, lca_squared, majority, lca_star



    ## LCA^2, Majority, and LCA* for each ORF




    exit()
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

    ## Build the LCA Star NCBI Tree
    print "Loading NCBI Tree..."
    lcastar = LCAStar(args["ncbi_tree"])
    print "Loaded."



    # set LCAStar parameters



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
        lineage = []
        for id in real_lineage:
            lineage.append(translate_to_prefered_name(id, ncbi_megan_map, lcastar))
        real_lineage = lineage
        line = "\t".join([contig, real, translate_to_prefered_name(lcastar.get_a_Valid_ID([taxon]),ncbi_megan_map, lcastar), "LCA_Star", str(lcastar.get_distance(taxon, real)), str(lcastar.wtd_distance(real, taxon)), ";".join(real_lineage[::-1])])
        output.write(line + "\n")
        taxon = lcastar.lca_majority(taxa_list)
        line = "\t".join([contig, real, translate_to_prefered_name(lcastar.get_a_Valid_ID([taxon]),ncbi_megan_map, lcastar), "Majority", str(lcastar.get_distance(taxon, real)), str(lcastar.wtd_distance(real, taxon)), ";".join(real_lineage[::-1])])
        output.write(line + "\n")
        taxon = lcastar.getTaxonomy(lca_list)
        line = "\t".join([contig, real, translate_to_prefered_name(lcastar.get_a_Valid_ID([taxon]),ncbi_megan_map, lcastar), "LCA_Squared", str(lcastar.get_distance(taxon, real)), str(lcastar.wtd_distance(real, taxon)), ";".join(real_lineage[::-1])])
        output.write(line + "\n")

    output.close()
# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
