#!/usr/bin/python
# File created on Nov 27 Jan 2012
from __future__ import division

__author__ = "Kishori M Konwar"
__copyright__ = "Copyright 2010, The metapaths Project"
__credits__ = ["r"]
__version__ = "1.0"
__maintainer__ = "Kishori M Konwar"
__status__ = "Release"

try:
     import pickle
     from cStringIO import StringIO
     import os
     from os import  makedirs, sys, remove, path
     import re
     from optparse import OptionParser, OptionGroup
     from LCAStar import *
     #from libs.python_modules.sysutil import getstatusoutput
except:
     print """ Could not load some user defined  module functions"""
     print """ """
     sys.exit(3)


objective = re.compile(r'Objective.*=\s*(\d*)')


usage=  sys.argv[0] + " --ncbi-tree-file ncbi_taxfile --taxalist taxa_list_file  -o output [--alpha majority --min_depth  min_depth_cutoff --min_reads min_read]"""


parser = OptionParser(usage)
parser.add_option("--ncbi-tree-file", dest="ncbi_tree_file", default = None,
                  help='name of the refseq names file [OPTIONAL]')

parser.add_option("--taxalist", dest="taxa_list", default = None,
                  help='list of taxons for contigs')

parser.add_option("--alpha", dest="alpha", type='float', default = 0.53, 
                  help='minimum fraction presetn to decide on a majority')

parser.add_option("--min_depth", dest="min_depth", type='int',  default = 3, 
                  help='minimum depth of the accepted taxon in the ncbi tree [OPTIONAL]')

parser.add_option("--min_reads", dest="min_reads", type = 'int', default = 5, 
                  help='minimum number of reads required after \"min_depth\" filtering [OPTIONAL]')

parser.add_option("-o", "--output", dest="output_file",
                  help='output file')



def create_annotation(namefile):

    file = 'ncbi_taxonomy_tree.txt'
    lca = SpeciesComputation(file)

    #taxonomy=lca.getTaxonomy(species)

def read_list_file_to_list(filename,list):
    namefile = open(filename, 'r')
    lines = namefile.readlines()
  
    for line in lines:
       line = line.strip()
       if line:
          list.append(line)


def read_taxonomy_list(filename, taxadict):
    try:
       namefile = open(filename, 'r')
    except:
       print "ERROR : Could not open file  " + filename
       sys.exit(1)
  
    line = True
    while line:
       line = namefile.readline()
       fields =  map(str.strip, line.split('\t') )
       if len(fields ) != 2:
          continue 
       if not fields[0] in taxadict:
           taxadict[fields[0]] = []

       taxadict[fields[0]].append(fields[1])

    try:
       namefile.close()
    except:
       print "ERROR : Could not close file  " + filename
       sys.exit(1)

   
def check_arguments(opts):
    if opts.ncbi_tree_file == None:
       print usage
       sys.exit(-1)

    if opts.taxa_list == None:
       print usage
       sys.exit(-1)


def main(argv): 

    (opts, args) = parser.parse_args()
    argv = check_arguments(opts)


    try:
       storagefileh = open('ncbi_tree_serialized', 'wb') 
    except:
       print 'ERROR : cannot open file to searialize the ncbi tree'
       sys.exit(-1)

    #p = pickle.Pickler(storagefileh)

    taxadict = {}
    read_taxonomy_list(opts.taxa_list, taxadict)

    lcastar = LCAStar(opts.ncbi_tree_file)
    lcastar.setLCAStarParameters(min_depth = opts.min_depth, alpha = opts.alpha, min_reads = opts.min_reads )
    #p.dump(lcastar)

    print 'done instatiating the lcastar'

    for contig in taxadict.keys():
       taxon = lcastar.lca_star(taxadict[contig])
       print contig + ' ' + taxon
       sys.exit(0)

    sys.exit(-1)



# the main function of metapaths
if __name__ == "__main__":
    main(sys.argv[1:])

