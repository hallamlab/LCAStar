#!/usr/bin/python

try:
     import sys, traceback
     import re
     import sys
#     from   metapaths_utils  import fprintf, printf, GffFileParser
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

    name_to_id={}
    id_to_name={}
    taxid_to_ptaxid = {}
    lca_min_score = 50
    lca_top_percent = 10 
    lca_min_support = 5
    results_dictionary = None
    def __init__(self, filename):
       taxonomy_file = open(filename, 'r')
       lines = taxonomy_file.readlines()
       taxonomy_file.close()

       for line in lines:
          if self.begin_pattern.search(line):
              continue
          fields =  [ x.strip()  for x in line.rstrip().split('\t')]
          if len(fields) !=3:
              continue
          self.name_to_id[str(fields[0])] = str(fields[1])
          self.id_to_name[str(fields[1])] = str(fields[0])
          self.taxid_to_ptaxid[str(fields[1])] = [ str(fields[2]), 0, 0]

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

    def translateNameToID(self, name):
       if not name in self.name_to_id:
           return None
       return self.name_to_id[name]

    def translateIdToName(self, id):
       if not id in self.id_to_name:
           return None
       return self.id_to_name[id]


    def getParentName(self, name):
       if not name in  self.name_to_id:  
          return None
       id = self.name_to_id[name]  
       pid = self.getParentTaxId(id)
       return self.translateIdToName( pid )


    def getParentTaxId(self, ID):
       if not ID in self.taxid_to_ptaxid:
          return None
       return self.taxid_to_ptaxid[ID][0]


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


    def getTaxonomy(self, name_groups):
         IDs = []
         for name_group in name_groups:
            id = self.get_a_Valid_ID(name_group)
            if id!=-1:
              IDs.append(id)
    
         consensus = self.get_lca(IDs)
         self.clear_cells(IDs)
         return consensus


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

    def set_results_dictionary(self, results_dictionary):
        self.results_dictionary= results_dictionary


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
 

