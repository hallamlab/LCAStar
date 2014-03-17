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

what_i_do = "Randomly select sequences to sample"
parser = argparse.ArgumentParser(description=what_i_do)

parser.add_argument('-i', dest='fasta_list', type=str, nargs='?',
             required=True, help='a list of the fasta files to sample from', default=None)                    
parser.add_argument('-o', dest='output_file', type=str, nargs='?',
             required=False, help='output file of simulated metagenome', default=None)
parser.add_argument('-l', dest='seq_length', type=int, nargs='?',
             required=False, help='length of sequence to sample', default=10000)
parser.add_argument('-s', dest='n_sequences_per_file', type=int, nargs='?',
             required=False, help='number of sample sequences', default=10)             
parser.add_argument('-n', dest='n_samples', type=int, nargs='?',
              required=False, help='number of sequences from list', default=100)

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
    
    # read in filenames
    sample_ids = [] # genomes to sample
    id_to_file = {} # ncbi_id to file map
    n_samples = args["n_samples"]
    
    file_list_h = open(args["fasta_list"], "r")
    lines = file_list_h.readlines() # assume not very big
    
    for l in lines:
        fields = l.split("\t")
        ncbi_id = fields[0].strip()
        ncbi_id = re.sub("\.[0-9]+$","", ncbi_id) #taking off version number
        f = fields[1].strip().strip("\n")
        id_to_file[ncbi_id] = f
        if ncbi_id not in sample_ids:
            sample_ids.append(ncbi_id)
    
    # randomize the sample_ids
    random.seed(12345) # set seed
    random.shuffle(sample_ids)
    
    # trunkate list to maximum size in n_samples
    if n_samples < len(sample_ids):
        sample_ids = sample_ids[0:n_samples]
    
    # open output file
    out_file = open(args["output_file"], "w")
    
    
    
    for i in sample_ids:
        reader = FastaReader(id_to_file[i])
        for record in reader:
            if len(record.sequence) > args["seq_length"]:
                for j in range(args["n_sequences_per_file"]):
                    # pick n_sequence_per_file random subsequences
                    ind = random.randint(0, len(record.sequence)-args["seq_length"])
                    end = ind + args["seq_length"]
                    name = "".join( [record.name , " ", "(", str(ind), "-", str(end),")"] )
                    out_file.write(name + "\n") # header with range sampled
                    out_file.write(record.sequence[ ind : end ] + "\n") # random sub-sequence
    
    out_file.close()

# the main function of metapaths
if __name__ == "__main__":
   main(sys.argv[1:])
