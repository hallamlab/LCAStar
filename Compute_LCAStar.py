#!python
"""
create_analysis_table.py

Created by Niels Hanson
Copyright (c) 2015 Steven J. Hallam Laboratory. All rights reserved.
"""

from __future__ import division

__author__ = "Niels W Hanson"
__copyright__ = "Copyright (c) 2015 Steven J. Hallam Laboratory. All rights reserved."
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
     import operator
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
    help='MetaPathways Output: Input mapping file in preprocessed directory (.mapping.txt).', default=None)
parser.add_argument('--ncbi_tree', dest='ncbi_tree', type=str, nargs='?', required=False,
    help='MetaPathways: NCBI tree file', default=None)
parser.add_argument('--ncbi_megan_map', dest='ncbi_megan_map', type=str, nargs='?', required=False,
    help='Preferred mapping for NCBI names (map.ncbi)', default=None)
parser.add_argument('-o', '--output', dest='output', type=str, nargs='?', required=False,
    help='output file of predicted taxonomies', default=None)
parser.add_argument('--orf_summary', dest='orf_summary', type=str, nargs='?',
    choices=['lca', 'besthit', 'orf_majority'], required=False, default='lca',
    help='ORF Summary method')
parser.add_argument('--contig_taxa_ref', dest='contig_taxa_ref', type=str, nargs='?', required=False,
    default=None, help='List of contig reference taxonomies (i.e., the known taxonomy)')
parser.add_argument('--sample_taxa_ref', dest='sample_taxa_ref', type=str, nargs='?', required=False,
    default=None, help='Name of the NCBI reference taxonomy. Hint: Put in double quotes')
parser.add_argument('--all_methods', dest='all_methods', action='store_true', required=False,
    default=None, help='Print all taxonomic estimation methods.')
parser.add_argument('--print_lineage', dest='print_lineage', action='store_true', required=False,
    default=None, help='Print full taxonomic linage instead of just final leaf taxonomy.')
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", required=False,
                    default=None, help="Verbose mode: prints all election counts for each contig.")

def translate_to_prefered_name(id, ncbi_megan_map, lcastar):
    # This maps an NCBI Taxonomy Database ID to the prefered MEGAN Taxonomy name and
    # reports the default name on the NCBI Taxonomy Database otherwise. 
    # id: id to be translated
    # ncbi_megan_map: dictionary of NCBI ID to MEGAN name
    # lcastar: instance of the LCAStar class
    id_str = str(id)
    if id_str in ncbi_megan_map:
        return ncbi_megan_map[id_str] + " (" + id_str + ")"
    else:
        res = lcastar.translateIdToName(id)
        if res:
            return res + " (" + id_str + ")"
        else:
            return "Unknown (" + id_str + ")"

def print_warning(message, error=False):
    if error:
        print "Error: " + message
    else:
        print "Warning: " + message


def clean_tab_lines(line):
    # Splits lines into fields by tabs, strips whitespace and end-of-line characters from each field.
    fields = line.split("\t")
    fields = map(str.strip, fields)
    fields = map(str.strip, fields, "\n")
    return fields

def create_contig2origin(mapping_filename):
    # dictionary/hash table mapping input sequences to their original header, e.g. >[Header]
    contig2origin = {}
    try:
        with open(mapping_filename, "r") as fh:
           for l in fh:
               fields = clean_tab_lines(l)
               read_name = fields[0]
               original_name = fields[1]
               contig2origin[read_name] = original_name
        return contig2origin
    except:
        print_warning("Could not open file.")

## global regular expression patterns
contig_pattern = re.compile("^(.*)_([0-9]*)$") # pulls out contig name from ORF annotations
taxonomy_pattern = re.compile("\[(.*)\]") # pulls out taxonomy between square brackets in taxonomic annotations

def verbose_writeout(args, contig_to_lca, contig_to_taxa_ref, sample_ref, lcastar, ncbi_megan_map):
    output_fh = None
    # if args['output']:
    #         output_fh = open(args['output'],"w")
    
    for contig in contig_to_lca:
        real = None
        
        # Prepare data structures
        orf_lcas = []
        simple_list = []
        taxa_to_orfs = {}
        real = None
        lca_squared_dist = None
        lca_squared_WTD = None
        majority_dist = None
        majority_WTD = None
        lca_star_dist = None
        lca_star_WTD = None
        
        for orf in contig_to_lca[contig]:
            lca = contig_to_lca[contig][orf]
            orf_lcas.append( [ lca ] )
            simple_list.append( lca )
            if lca not in taxa_to_orfs:
                taxa_to_orfs[lca] = []
            taxa_to_orfs[lca].append(orf)
        
        # Calculate statistics and p-values
        lca_squared_id = lcastar.getTaxonomy( orf_lcas, return_id=True )
        majority = lcastar.simple_majority( simple_list) # for p-val
        majority_id = lcastar.simple_majority( simple_list, return_id=True )
        majority_p = lcastar.calculate_pvalue(simple_list, majority)
        lca_star_id, lca_star_p = lcastar.lca_star( simple_list, return_id=True )
        
        # Calcualte distances frome expected if needed
        if contig_to_taxa_ref or sample_ref:
            lca_squared = lcastar.getTaxonomy( orf_lcas )
            lca_star, lca_star_p = lcastar.lca_star( simple_list )
            
            if contig_to_taxa_ref:
                real = contig_to_taxa_ref[contig]
            else:
                real = sample_ref
            lca_squared_dist = str(lcastar.get_distance(lca_squared, real))
            lca_squared_WTD = str(lcastar.wtd_distance(lca_squared, real))
            majority_dist = str(lcastar.get_distance(majority, real))
            majority_WTD = str(lcastar.wtd_distance(majority, real))
            lca_star_dist = str(lcastar.get_distance(lca_star, real))
            lca_star_WTD = str(lcastar.wtd_distance(lca_star, real))
        
        # Get lineages 
        # TODO: make list only have individual taxa if not printing full lineages
        lca_squared_lineage_ids = lcastar.get_lineage(lca_squared_id)
        majority_lineage_ids = lcastar.get_lineage(majority_id)
        lca_star_lineage_ids = lcastar.get_lineage(lca_star_id)
        
        # Translate ids to preferred names and join by ';'
        lca_squared_lineage =  ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lca_squared_lineage_ids[::-1]])
        majority_lineage = ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in majority_lineage_ids[::-1]])
        lca_star_lineage = ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lca_star_lineage_ids[::-1]])
        
        print "Contig:", contig,
        if real:
            print "Origin:", real
        
        print ""
        for taxa in taxa_to_orfs:
            orf_ids = []
            for x in taxa_to_orfs[taxa]:
                orf_ids.append(contig + "_" + x)
            print taxa, ":", len(taxa_to_orfs[taxa]), orf_ids
        
        print ""
        print "LCASquared:", lca_squared_lineage
        if contig_to_taxa_ref or sample_ref:
            print "Distances:(" + str(lca_squared_dist) + ", " + str(lca_squared_WTD) + ")"
        print "Majority:", majority_lineage, "p =", majority_p
        if contig_to_taxa_ref or sample_ref:
            print "Distances:(" + str(majority_dist) + ", " + str(majority_WTD) + ")"
        print "LCAStar:", lca_star_lineage, "p =", lca_star_p
        if contig_to_taxa_ref or sample_ref:
            print "Distances:(" + str(lca_star_dist) + ", " + str(lca_star_WTD) + ")"
        print "--------"
    exit()
            # 
            # 
            #     ## calculate distances
            #     real = None
            #     lca_squared_dist = None
            #     lca_squared_WTD = None
            #     majority_dist = None
            #     majority_WTD = None
            #     lca_star_dist = None
            #     lca_star_WTD = None
            # 
            #     if contig_to_taxa_ref or sample_ref:
            #         if contig_to_taxa_ref:
            #             real = contig_to_taxa_ref[contig]
            #         else:
            #             real = sample_ref
            #         lca_squared_dist = str(lcastar.get_distance(lca_squared, real))
            #         lca_squared_WTD = str(lcastar.wtd_distance(lca_squared, real))
            #         majority_dist = str(lcastar.get_distance(majority, real))
            #         majority_WTD = str(lcastar.wtd_distance(majority, real))
            #         lca_star_dist = str(lcastar.get_distance(lca_star, real))
            #         lca_star_WTD = str(lcastar.wtd_distance(lca_star, real))
            #     
            #     # construct line based on cases
            #     if all_methods:
            #         if contig_to_taxa_ref or sample_ref:
            #             line = "\t".join(map(str, [contig, lca_star, lca_star_p, lca_star_dist, lca_star_WTD,
            #                                  majority, majority_p, majority_dist, majority_WTD,
            #                                  lca_squared, lca_squared_dist, lca_squared_WTD, real]))
            #         else:
            #             line = "\t".join(map(str, [contig, lca_star, lca_star_p,
            #                                        majority, majority_p, lca_squared]))
            #     else:
            #         if contig_to_taxa_ref or sample_ref:
            #             line = "\t".join([contig, lca_star, lca_star_p, lca_star_dist, lca_star_WTD, real ])
            #         else:
            #             line = "\t".join([contig, lca_star, lca_star_p])
            # 
            #     # print out line
            #     if output_fh:
            #         output_fh.write(line + "\n")
            #     else:
            #         print line
            # if output_fh:
            #     output_fh.close()
    

def main(argv):
    # parse arguments
    args = vars(parser.parse_args())

    ## read input mapping file (.mapping.txt) to create read2origin map
    contig2origin = {} # dictionary/hash table mapping input sequences to their original header, e.g. >[Header]
    contig2origin = create_contig2origin(args["mapping_file"])

    ## create preferred mapping
    ncbi_megan_map = {} # hash map from given taxonomy to preferred one used by megan
    with open(args["ncbi_megan_map"], 'r') as meganfile:
        for line in meganfile:
             fields = line.split("\t")
             fields = map(str.strip, fields)
             ncbi_megan_map[fields[0]] = fields[1]

    # contig_taxa_data structure
    contig_to_taxa = {}

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
                        bitscore = fields[3]
                        contig_to_taxa[contig][orf].append( (taxa, int(bitscore)) )
                    else:
                        continue
                else:
                    continue

    ## load contig references (if available or applicable)

    # read contig taxa reference if available
    contig_to_taxa_ref = None
    if args["contig_taxa_ref"]:
        contig_to_taxa_ref = {}
        if args["contig_taxa_ref"]:
            with open(args["contig_taxa_ref"], "r") as fh:
                for l in fh:
                    fields = clean_tab_lines(l)
                    contig_id = fields[0]
                    contig_origin = fields[1]
                    contig_to_taxa_ref[contig_id] = contig_origin

    # all contigs hypothetically have the same reference origin (i.e., single cells)
    sample_ref = None
    if args["sample_taxa_ref"]:
        sample_ref = args["sample_taxa_ref"]

    ## Build the LCA Star NCBI Tree
    lcastar = LCAStar(args["ncbi_tree"])
    lcastar.setLCAStarParameters(min_depth = 1, alpha = 0.51, min_reads = 1)

    ## Calculate LCA for each ORF
    contig_to_lca = {}
    for contig in contig_to_taxa:
        for orf in contig_to_taxa[contig]:
            if contig not in contig_to_lca:
                contig_to_lca[contig] = {}
            if orf not in contig_to_lca[contig]:
                contig_to_lca[contig][orf] = None
            contig_taxas = contig_to_taxa[contig][orf]
            if len(contig_taxas) == 0:
                contig_to_lca[contig][orf] = "root"
            else:
                if args['orf_summary'] == 'besthit':
                    contig_taxas.sort(key=operator.itemgetter(1), reverse=True)
                    best_blast_taxa = contig_taxas[0][0]
                    contig_to_lca[contig][orf] = best_blast_taxa
                elif args['orf_summary'] == 'orf_majority':
                    majority_list = []
                    for t in contig_taxas:
                        majority_list.append(t[0])
                    #TODO Update to check for alterantive taxonomy names
                    contig_to_lca[contig][orf] = lcastar.simple_majority(majority_list)
                else:
                    lca_list = [] # create a list of lists for LCA calculation
                    for t in contig_taxas:
                        lca_list.append([t[0]])
                    contig_to_lca[contig][orf] = lcastar.getTaxonomy(lca_list)

    ## calculate taxonomy statistics LCA,  for each ORF
    # contig_to_taxa = {}

    ## LCA^2, Majority, and LCA* for each ORF
    if args['all_methods']:
        if contig_to_taxa_ref or sample_ref:
            header = "\t".join(["Contig", "LCAStar", "LCAStar_p", "LCAStar_dist", "LCAStar_WTD",
                         "Majority", "Majority_p", "Majority_dist", "Majority_WTD",
                         "LCASquared", "LCASquared_dist", "LCASquared_WTD", "Original"])
        else:
            header = "\t".join(["Contig", "LCAStar", "LCAStar_p",
                         "Majority", "Majority_p", "LCASquared"])
    else:
        if contig_to_taxa_ref or sample_ref:
            header = "\t".join(["Contig", "LCAStar", "LCAStar_p", "LCAStar_dist", "LCAStar_WTD", "Original" ])
        else:
            header = "\t".join(["Contig", "LCAStar", "LCAStar_p"])

    output_fh = None
    if args['output']:
        output_fh = open(args['output'],"w")
        output_fh.write(header + "\n")
    else:
        print header
    
    verbose_writeout(args, contig_to_lca, contig_to_taxa_ref, sample_ref, lcastar, ncbi_megan_map)
    
    for contig in contig_to_lca:
        orf_lcas = []
        simple_list = []
        for orf in contig_to_lca[contig]:
            lca = contig_to_lca[contig][orf]
            orf_lcas.append( [ lca ] )
            simple_list.append( lca )

        ## calculate statistics and p-values
        lca_squared = lcastar.getTaxonomy( orf_lcas )
        majority = lcastar.simple_majority( simple_list )
        majority_p = lcastar.calculate_pvalue(simple_list, majority)
        lca_star, lca_star_p = lcastar.lca_star( simple_list )

        ## calculate distances
        real = None
        lca_squared_dist = None
        lca_squared_WTD = None
        majority_dist = None
        majority_WTD = None
        lca_star_dist = None
        lca_star_WTD = None

        if contig_to_taxa_ref or sample_ref:
            if contig_to_taxa_ref:
                real = contig_to_taxa_ref[contig]
            else:
                real = sample_ref
            lca_squared_dist = str(lcastar.get_distance(lca_squared, real))
            lca_squared_WTD = str(lcastar.wtd_distance(lca_squared, real))
            majority_dist = str(lcastar.get_distance(majority, real))
            majority_WTD = str(lcastar.wtd_distance(majority, real))
            lca_star_dist = str(lcastar.get_distance(lca_star, real))
            lca_star_WTD = str(lcastar.wtd_distance(lca_star, real))
        
        # construct line based on cases
        if args['all_methods']:
            if contig_to_taxa_ref or sample_ref:
                line = "\t".join(map(str, [contig, lca_star, lca_star_p, lca_star_dist, lca_star_WTD,
                                     majority, majority_p, majority_dist, majority_WTD,
                                     lca_squared, lca_squared_dist, lca_squared_WTD, real]))
            else:
                line = "\t".join(map(str, [contig, lca_star, lca_star_p,
                                           majority, majority_p, lca_squared]))
        else:
            if contig_to_taxa_ref or sample_ref:
                line = "\t".join([contig, lca_star, lca_star_p, lca_star_dist, lca_star_WTD, real ])
            else:
                line = "\t".join([contig, lca_star, lca_star_p])

        # print out line
        if output_fh:
            output_fh.write(line + "\n")
        else:
            print line
    if output_fh:
        output_fh.close()


    exit()






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
