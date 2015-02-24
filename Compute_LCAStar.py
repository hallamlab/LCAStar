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
parser.add_argument('-a', '--all_methods', dest='all_methods', action='store_true', required=False,
    default=None, help='Print all taxonomic estimation methods.')
parser.add_argument('-l','--print_lineage', dest='print_lineage', action='store_true', required=False,
    default=None, help='Print full taxonomic linage instead of just final leaf taxonomy.')
parser.add_argument("-v", "--verbose", dest="verbose", action="store_true", required=False,
                    default=None, help="Verbose mode: prints all election counts for each contig.")
parser.add_argument("--alpha", dest="alpha", type=float, required=False,
                    default=0.51, help="Alpha-majority threshold.")

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
        res = lcastar.translateIdToName(id_str)
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

def printline(line, output_fh=None):
    # print or write line depending on presence of output_fh
    if output_fh:
        output_fh.write(line + "\n")
    else:
        print line

def writeout(args, contig_to_lca, contig_to_taxa_ref, sample_ref, lcastar, ncbi_megan_map):
    output_fh = None
    if args['output']:
       output_fh = open(args['output'],"w")
    
    if args["verbose"]:
        for contig in contig_to_lca:
            # Prepare data structures
            real = None
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
            lca_squared_lineage_ids = None
            majority_lineage_ids = None
            lca_star_lineage_ids = None
            
            # Collect taxonomy for each ORF
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
        
            # Calcualte distances from expected if needed
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
            lca_squared_lineage_ids = lcastar.get_lineage(lca_squared_id)
            majority_lineage_ids = lcastar.get_lineage(majority_id)
            lca_star_lineage_ids = lcastar.get_lineage(lca_star_id)
            
            # Remove all but last member if full-lineage not requested
            if not args["print_lineage"]:
                lca_squared_lineage_ids = [lca_squared_lineage_ids[0]]
                majority_lineage_ids = [majority_lineage_ids[0]]
                lca_star_lineage_ids = [lca_star_lineage_ids[0]]
        
            # Translate ids to preferred names and join by ';'
            lca_squared_lineage =  ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lca_squared_lineage_ids[::-1]])
            majority_lineage = ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in majority_lineage_ids[::-1]])
            lca_star_lineage = ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lca_star_lineage_ids[::-1]])
            
            # Print verbose output for contig
            line = "Contig: " + str(contig)
            printline(line, output_fh)
                
            if real:
                line = "Origin: " + str(real)
                printline(line, output_fh)
            
            printline("", output_fh)
             
            for taxa in taxa_to_orfs:
                orf_ids = []
                for x in taxa_to_orfs[taxa]:
                    orf_ids.append(contig + "_" + x)
                taxa_alt = translate_to_prefered_name(lcastar.translateNameToID(taxa), ncbi_megan_map, lcastar)
                line = "".join([taxa_alt, ": ", str(len(taxa_to_orfs[taxa])), " [",",".join(orf_ids),"]"])
                printline(line, output_fh)
        
            printline("", output_fh)
            
            # Construct line based on cases
            if args['all_methods']:
                line = "".join(["LCASquared: ", lca_squared_lineage])
                printline(line, output_fh)
                if contig_to_taxa_ref or sample_ref:
                    line = "Distances: (" + str(lca_squared_dist) + ", " + str(lca_squared_WTD) + ")"
                    printline(line, output_fh)
                line = "".join(["Majority: ", majority_lineage, " p=", str(majority_p)])
                printline(line, output_fh)
                if contig_to_taxa_ref or sample_ref:
                    line = "Distances:(" + str(majority_dist) + ", " + str(majority_WTD) + ")"
                    printline(line, output_fh)
                line = "".join(["LCAStar: ", lca_star_lineage, " p=", str(lca_star_p)])
                printline(line, output_fh)
                if contig_to_taxa_ref or sample_ref:
                    line = "Distances: (" + str(lca_star_dist) + ", " + str(lca_star_WTD) + ")"
                    printline(line, output_fh)
            else:
                line = "".join(["LCAStar: ", lca_star_lineage, " p=", str(lca_star_p)])
                printline(line, output_fh)
                if contig_to_taxa_ref or sample_ref:
                    line = "Distances: (" + str(lca_star_dist) + ", " + str(lca_star_WTD) + ")"
                    printline(line, output_fh)
            printline("--------", output_fh)
    else:
        # print table header
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
        
        printline(header, output_fh)
        
        for contig in contig_to_lca:
            # Prepare data structures
            real = None
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
            lca_squared_lineage_ids = None
            majority_lineage_ids = None
            lca_star_lineage_ids = None
            
            orf_lcas = []
            simple_list = []
            for orf in contig_to_lca[contig]:
                lca = contig_to_lca[contig][orf]
                orf_lcas.append( [ lca ] )
                simple_list.append( lca )

            # Calculate statistics and p-values
            lca_squared_id = lcastar.getTaxonomy( orf_lcas, return_id=True )
            majority = lcastar.simple_majority( simple_list) # for p-val
            majority_id = lcastar.simple_majority( simple_list, return_id=True )
            majority_p = lcastar.calculate_pvalue(simple_list, majority)
            lca_star_id, lca_star_p = lcastar.lca_star( simple_list, return_id=True )

            # Calcualte distances from expected if needed
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
            lca_squared_lineage_ids = lcastar.get_lineage(lca_squared_id)
            majority_lineage_ids = lcastar.get_lineage(majority_id)
            lca_star_lineage_ids = lcastar.get_lineage(lca_star_id)
            
            # Remove all but last member if full-lineage not requested
            if args["print_lineage"] == None:
                lca_squared_lineage_ids = [lca_squared_lineage_ids[0]]
                majority_lineage_ids = [majority_lineage_ids[0]]
                lca_star_lineage_ids = [lca_star_lineage_ids[0]]
            
            # Translate ids to preferred names and join by ';'
            lca_squared_lineage =  ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lca_squared_lineage_ids[::-1]])
            majority_lineage = ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in majority_lineage_ids[::-1]])
            lca_star_lineage = ";".join([translate_to_prefered_name(x, ncbi_megan_map, lcastar) for x in lca_star_lineage_ids[::-1]])
        
            # Construct line based on cases
            if args['all_methods']:
                if contig_to_taxa_ref or sample_ref:
                    line = "\t".join(map(str, [contig, lca_star_lineage, lca_star_p, lca_star_dist, lca_star_WTD,
                                         majority_lineage, majority_p, majority_dist, majority_WTD,
                                         lca_squared_lineage, lca_squared_dist, lca_squared_WTD, real]))
                else:
                    line = "\t".join(map(str, [contig, lca_star_lineage, lca_star_p,
                                               majority_lineage, majority_p, lca_squared_lineage]))
            else:
                if contig_to_taxa_ref or sample_ref:
                    line = "\t".join([contig, lca_star_lineage, lca_star_p, lca_star_dist, lca_star_WTD, real ])
                else:
                    line = "\t".join([contig, lca_star_lineage, lca_star_p])

            # Print out line
            if output_fh:
                output_fh.write(line + "\n")
            else:
                print line
    if output_fh:
        output_fh.close()
    
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

    # Read blast table
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

    ## Load contig references (if available or applicable)

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
    lcastar.setLCAStarParameters(min_depth = 1, alpha = args["alpha"], min_reads = 1)

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
    writeout(args, contig_to_lca, contig_to_taxa_ref, sample_ref, lcastar, ncbi_megan_map)

# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
