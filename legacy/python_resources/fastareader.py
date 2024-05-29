
import re

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

