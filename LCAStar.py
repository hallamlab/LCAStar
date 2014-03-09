#!/usr/bin/python

try:
     import sys, traceback
     import re
     import sys
     from math import log
except:
     print """ Could not load some user defined  module functions"""
     print """ Make sure your typed \"source MetaPathwaysrc\""""
     print """ """
     print traceback.print_exc(10)
     sys.exit(3)



def copyList(a, b): 
    [ b.append(x) for x in a ] 

class LCAStar:
    begin_pattern = re.compile("#")

    # a readable taxon name to numeric string id map as ncbi
    name_to_id={}

    # a readable taxon ncbi tax id to name map
    id_to_name={}

    # this is the tree structure in a id to parent map, you can traverse it to go to the root
    taxid_to_ptaxid = {}

    # this is the tree structure in a parent to child id map, you can use it to traverse the tree downwards
    # ptaxid_to_taxid[ptaxid] = [ cid1, cid2, ...cidk]
    ptaxid_to_taxid = {}

    # a map from id to value, which has the S = sum n,  value for each id
    id_to_R={}
    # a map from id to value, which has the S = sum n,  value for each id
    id_to_S={}
    # a map from id to value, which has the L = sum n log n,  value for each id
    id_to_L={}
    # a map from id to value, which has the entropy H value for each id
    id_to_H={}

    # a map to keep track of visited nodes
    id_to_V={}


    lca_min_score = 50   # an LCA parameter for min score for a hit to be considered
    lca_top_percent = 10    # an LCA param to confine the hits to within the top hits score upto the top_percent% 
    lca_min_support = 5   # a minimum number of reads in the sample to consider a taxon to be present

    lca_star_min_reads = 10
    lca_star_min_depth = 3
    lca_star_alpha = 0.53

    results_dictionary = None

    # initialize with the ncbi tree file 
    def __init__(self, filename):
       taxonomy_file = open(filename, 'r')
       lines = taxonomy_file.readlines()
       taxonomy_file.close()

       for line in lines:
          if self.begin_pattern.search(line):
              continue
          fields =  [ str(x.strip())  for x in line.rstrip().split('\t')]
          if len(fields) !=3:
              continue
          self.name_to_id[fields[0]] = fields[1]
          self.id_to_name[fields[1]] = fields[0]
          # the taxid to ptax map has for each taxid a corresponding 3-tuple
          # the first location is the pid, the second is used as a counter for 
          # lca while a search is traversed up the tree and the third is used for
          # the min support
          self.taxid_to_ptaxid[fields[1]] = [fields[2], 0, 0]

          if not fields[2] in self.ptaxid_to_taxid:
              self.ptaxid_to_taxid[fields[2]]={}

          if not (fields[2]=='1' and fields[1]=='1'):
             if fields[1]==fields[2]:
                 print "ERROR what did you just do!. You should  have never  gotten here!"
                 sys.exit(0)

             self.ptaxid_to_taxid[fields[2]][fields[1]] = False


    def setParameters(self, min_score, top_percent, min_support):
       self.lca_min_score = min_score
       self.lca_top_percent =top_percent
       self.lca_min_support = min_support
         
    def sizeTaxnames(self ):
         return len(self.name_to_id)


    def sizeTaxids(self):
         return len(self.taxid_to_ptaxid)
          
    def get_a_Valid_ID(self, name_group):
        for name in name_group:
           if name in self.name_to_id:
               return  self.name_to_id[name]
        return -1

    # given a taxon name it returns the correcponding unique ncbi tax id
    def translateNameToID(self, name):
       if not name in self.name_to_id:
           return None
       return self.name_to_id[name]

    # given a taxon id to taxon name map
    def translateIdToName(self, id):
       if not id in self.id_to_name:
           return None
       return self.id_to_name[id]


    # given a name it returns the parents name
    def getParentName(self, name):
       if not name in  self.name_to_id:  
          return None
       id = self.name_to_id[name]  
       pid = self.getParentTaxId(id)
       return self.translateIdToName( pid )


    # given a ncbi tax id returns the parents tax id
    def getParentTaxId(self, ID):
       if not ID in self.taxid_to_ptaxid:
          return None
       return self.taxid_to_ptaxid[ID][0]


    # given a set of ids it returns the lowest common ancenstor 
    # without caring about min support
    # here LCA for a set of ids are computed as follows
    # first we consider one ID at a time
    #   for each id we traverse up the ncbi tree using the id to parent id map
    #   at the same time increasing the count on the second value of the 3-tuple 
    #   note that at the node where all the of the individual ids ( limit in number)
    #   converges the counter matches the limit for the first time, while climbing up. 
    #   This also this enables us to  make the selection of id arbitrary 
    def get_lca(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               self.taxid_to_ptaxid[tid][1]+=1
               if self.taxid_to_ptaxid[tid][1]==limit:
                  return  self.id_to_name[tid]  
               tid = self.taxid_to_ptaxid[tid][0]

        return "root"

    def update_taxon_support_count(self, taxonomy):
         id = self.get_a_Valid_ID( [taxonomy ])
         tid = id 
         while( tid in self.taxid_to_ptaxid and tid !='1' ):
               self.taxid_to_ptaxid[tid][2]+=1
               tid = self.taxid_to_ptaxid[tid][0]

    def get_supported_taxon(self, taxonomy):
         id = self.get_a_Valid_ID( [taxonomy ])
         tid = id 
         while( tid in self.taxid_to_ptaxid and tid !='1' ):
            if self.lca_min_support > self.taxid_to_ptaxid[tid][2] :
                tid = self.taxid_to_ptaxid[tid][0]
            else:
                return self.translateIdToName(tid)

         return  self.translateIdToName(tid)
    
    # need to call this to clear the counts of reads at every node      
    def clear_cells(self, IDs):
        limit = len(IDs)
        for id in IDs:
           tid = id 
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
               #if self.taxid_to_ptaxid[tid][1]==0:
               #   return  self.id_to_name[tid]  
               self.taxid_to_ptaxid[tid][1]=0
               tid = self.taxid_to_ptaxid[tid][0]
        return ""


    #given a set of sets of names it computes an lca 
    # in the format [ [name1, name2], [name3, name4,....namex] ...]
    # here name1 and name2 are synonyms and so are name3 through namex
    def getTaxonomy(self, name_groups):
         IDs = []
         for name_group in name_groups:
            id = self.get_a_Valid_ID(name_group)
            if id!=-1:
              IDs.append(id)
    
         consensus = self.get_lca(IDs)
         self.clear_cells(IDs)
         return consensus


    # extracts taxon names for a refseq annotation
    def get_species(self, hit):
       if not 'product' in hit: 
           return None
       species = []
       try:
           m = re.findall(r'\[([^\[]+)\]', hit['product'])
           if m != None:
             copyList(m,species)
       except:
             return None
   
       if species:
          return species
       else:
          return None
 
    # used for optimization
    def set_results_dictionary(self, results_dictionary):
        self.results_dictionary= results_dictionary

    # this returns the megan taxonomy, i.e., it computes the lca but at the same time
    # takes into consideration the parameters, min score, min support and top percent
    def getMeganTaxonomy(self, orfid):
         #compute the top hit wrt score
         names = []
         species = []
         if 'refseq' in self.results_dictionary:
            if orfid in self.results_dictionary['refseq']:
                 
               top_score = 0 
               for hit in self.results_dictionary['refseq'][orfid]:
                  if hit['bitscore'] >= self.lca_min_score and hit['bitscore'] >= top_score:
                     top_score = hit['bitscore']

               for hit in self.results_dictionary['refseq'][orfid]:
                  if (100-self.lca_top_percent)*top_score/100 < hit['bitscore']:
                     names = self.get_species(hit)
                     #if 'MD_2_95' == orfid:
                     #  for hit in self.results_dictionary['refseq'][orfid]:
                     #     print  orfid  + ':' + str(names)
                     #  else:
                     #     print orfid  + ':' + str([])
                     if names:
                       species.append(names) 

         taxonomy = self.getTaxonomy(species)
         meganTaxonomy = self.get_supported_taxon( taxonomy)
         return meganTaxonomy
 

    # this is use to compute the min support for each taxon in the tree
    # this is called before the  getMeganTaxonomy
    def compute_min_support_tree(self, annotate_gff_file, pickorfs):
        gffreader = GffFileParser(annotate_gff_file)
        try:
           for contig in  gffreader:
              for orf in  gffreader.orf_dictionary[contig]:
                 if not orf['id'] in pickorfs:
                     continue
                 taxonomy = None
                 species = []
                 if 'refseq' in self.results_dictionary:
                   if orf['id'] in self.results_dictionary['refseq']:
                       #compute the top hit wrt score
                       top_score = 0 
                       for hit in self.results_dictionary['refseq'][orf['id']]:
                          if hit['bitscore'] >= self.lca_min_score and hit['bitscore'] >= top_score:
                            top_score = hit['bitscore']
       
                       for hit in self.results_dictionary['refseq'][orf['id']]:
                          if (100-self.lca_top_percent)*top_score/100 < hit['bitscore']:
                             names = self.get_species(hit)
                             if names:
                               species.append(names) 
                 taxonomy=self.getTaxonomy(species)
                 self.update_taxon_support_count(taxonomy)
                 pickorfs[orf['id']] = taxonomy
        except:
           print "ERROR : Cannot read annotated gff file "
          

    
    def taxon_depth(self, taxon):
        id = self.translateNameToID(taxon)
        if id==None:
           return 0

        tid = id 
        depth = 0
        #climb up the tree from the taxon to the root 
        # the number of climbing steps is the depth
        while( tid in self.taxid_to_ptaxid and tid !='1' ):
             tid = self.taxid_to_ptaxid[tid][0]
             depth +=1

        return depth

    def filter_taxa_list(self, taxalist):
        #filter baed on depth 
        newlist = []
        for taxon in taxalist:
            depth  = self.taxon_depth(taxon)
            if depth < self.lca_star_min_depth:
               continue
            newlist.append(taxon)

        #filter based on min_reads /decide if we should return 'root'
        # or compute the taxon using lca_star
        if len(newlist) < self.lca_star_min_reads:
            return None 
        else:
            return newlist

    def setLCAStarParameters(self,  min_depth = 3, alpha = 0.53,  min_reads = 10 ):
        self.lca_star_min_reads = min_reads
        self.lca_star_min_depth = min_depth
        self.lca_star_alpha = alpha
        
    def __read_counts(self, taxalist):
        read_counts = {}
        Total = 0
        for taxon in taxalist:
           if not taxon in read_counts:
             read_counts[taxon] = 0
           read_counts[taxon] += 1
           Total +=1
        return read_counts, Total
         
        # find the taxon with the highest count but also has count higher than the 
        # majority threshold
 
    def lca_majority(self, taxalist):
        taxalist = self.filter_taxa_list(taxalist)

        majority = 'all'
        if taxalist==None:
           return majority

        majority = self.__lca_majority(taxalist)

        if majority==None:
           return 'all'

        return majority
        

    def __lca_majority(self, taxalist):

        # create the read counts
        read_counts, Total = self.__read_counts(taxalist)
         
        # find the taxon with the highest count but also has count higher than the 
        # majority threshold
        majority = None
        maxcount =0
        for taxon in taxalist: 
           if maxcount < read_counts[taxon] and Total*self.lca_star_alpha < read_counts[taxon]:
               maxcount = read_counts[taxon]
               majority = taxon

        # majority exists 
        if majority!= None :
          return majority

        return None

    def __color_tree(self, read_counts):
        for taxon, value in read_counts.iteritems():
           id = self.translateNameToID(taxon)
           tid = id 
           #climb up the tree from the taxon to the root 
           #  and mark the parent to child structure with True
           while( tid in self.taxid_to_ptaxid and tid !='1' ):
              pid = self.taxid_to_ptaxid[tid][0]
              self.ptaxid_to_taxid[pid][tid] = True
              tid = pid

    def __annotate_tree_counts(self, read_counts):
        for taxon, value in read_counts.iteritems():
           id = self.translateNameToID(taxon)
           if id==None:
               continue
           self.id_to_R[id] = value


    def __decolor_tree(self):
        S = ['1']
        while len(S) >0:
            id = S.pop()

            C = []
            if id in self.ptaxid_to_taxid:
              C = self.ptaxid_to_taxid[id].keys()
           
            for child in C:
               if self.ptaxid_to_taxid[id][child]:
                  self.ptaxid_to_taxid[id][child] = False
                  S.append(child)

           
    def __create_majority(self, root, read_name_counts):
        read_counts = {}
        Total = 0
        for  taxon, count in read_name_counts.iteritems():
            id = self.translateNameToID(taxon)
            read_counts[id] = count
            Total += count

        candidate = ['1', 10000000.00 ]

        Stack = [root]
        while len(Stack) >0:
            id = Stack.pop()
            # calculate here
            if id in self.id_to_V:
               C = []
               if id in self.ptaxid_to_taxid:
                  C = self.ptaxid_to_taxid[id].keys()
               # I am coming up
               self.id_to_H[id] = 0

               if id in read_counts:
                  self.id_to_S[id] = float(read_counts[id])
                  self.id_to_L[id] = float(read_counts[id])*log(float(read_counts[id]))
               else:
                  self.id_to_S[id] = 0 
                  self.id_to_L[id] = 0 

               for child in C:
                 if self.ptaxid_to_taxid[id][child]:
                   self.id_to_S[id] += self.id_to_S[child] 
                   self.id_to_L[id] += self.id_to_L[child] 

               self.id_to_H[id] = -(self.id_to_L[id]/self.id_to_S[id] - log(self.id_to_S[id]))

               if self.id_to_S[id] > Total*self.lca_star_alpha:
                   if candidate[1] > self.id_to_H[id]:
                      candidate[0] = id
                      candidate[1] = self.id_to_H[id]

               #compute entropy calculation with S[id]

            else: # going down
               self.id_to_V[id] = True
               C = []
               if id in self.ptaxid_to_taxid:
                  C = self.ptaxid_to_taxid[id].keys()

               Stack.append(id)
               for child in C:
                 if self.ptaxid_to_taxid[id][child]:
                    Stack.append(child)
        return candidate[0] 

    def __clear_lca_star_data_structure(self):
        self.id_to_R={}
        self.id_to_S={}
        self.id_to_L={}
        self.id_to_H={}
        self.id_to_V={}

        
    def lca_star(self, taxalist):
        taxalist = self.filter_taxa_list(taxalist)

        if taxalist==None:
           return 'all'
        
        majority = self.__lca_majority(taxalist) 

        if majority != None:
           return majority 

        read_counts, Total = self.__read_counts(taxalist)
#        for key, value in read_counts.iteritems():
#           print key + ' ' + str(value) 

        self.__annotate_tree_counts(read_counts)
        self.__color_tree(read_counts)
        resultid = self.__create_majority('1', read_counts) 
        self.__clear_lca_star_data_structure()
        resulttaxon = self.translateIdToName( str(resultid ))
        self.__decolor_tree()
        return resulttaxon
