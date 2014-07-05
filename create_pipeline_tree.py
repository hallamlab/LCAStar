#!python
"""
create_pipeline_tree.py

Created by Niels Hanson
Copyright (c) 2013 Steven J. Hallam Laboratory. All rights reserved.
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

what_i_do = "Creates an NCBI Tree Compatible with MetaPathways."
parser = argparse.ArgumentParser(description=what_i_do)

parser.add_argument('--ncbi_names', dest='ncbi_names', type=str, nargs='?',
          required=True, help='NCBI Taxonomy tree names.dmp file', default=None)
parser.add_argument('--ncbi_nodes', dest='ncbi_nodes', type=str, nargs='?',
          required=True, help='NCBI Taxonomy tree nodes.dmp file', default=None)                
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
          required=False, help='name out output tree file', default=None)

def main(argv):
   args = vars(parser.parse_args())
   
   names = {}
   fh = open(args['ncbi_names'], "r")
   l = fh.readline()
   while l:
      fields = l.split("|")
      fields = map(str.strip, fields, "\t")
      fields = map(str.strip, fields)
      if fields[0] not in names:
          names[fields[0]] = []
      names[fields[0]].append(fields[1])
      l = fh.readline()
   fh.close()
   
   fh = open(args['ncbi_nodes'], "r")
   l = fh.readline()
   output = open(args['output_file'], "w")
   while l:
       fields = l.split("|")
       fields = map(str.strip, fields, "\t")
       child_id = fields[0]
       parent_id = fields[1]
       if child_id in names:
           while names[child_id]:
               output.write("\t".join([names[child_id].pop(), child_id, parent_id]) + "\n")
       else:
           print "Error: missing child"
       l = fh.readline()
   fh.close()
   output.close()
   


# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
