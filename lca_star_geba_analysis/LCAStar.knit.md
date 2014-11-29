---
title: "LCA Star"
output:
   html_document: 
      toc: true
      theme: readable
      highlight: default
---

**Goal**: To demonstrate the effectiveness of the current implementation of LCA* for taxonomic assignment with respect to existing methods of $LCA^2$, Contig LCA, and Simple Voting procedures with respect to classification error and sum of squares of NCBI Tree Distance from prediction and actual taxonomy.

## Simulation using MetaSim and NCBI Genomes

* Downloaded `all.fna.tar.gz` and `summary.txt` from the NCBI's ftp server:ftp://ftp.ncbi.nlm.nih.gov/genomes/ Mar 12, 2014
* `summary.txt` contains some metadata for a number of these genomes:

```
Accession  GenbankAcc	Length	Taxid	ProjectID	TaxName	Replicon	Create Date	Update_Date
NC_000907.1	L42023.1	1830138	71421	57771	Haemophilus influenzae Rd KW20	chromosome 	Oct 19 2001	Sep 11 2013  4:30:20:076PM
```

* this kind of information would be useful for an analysis so I cross referenced all the Accession IDs with the actual `.fna` files that I had in `all.fna.tar.gz` with the following shell command:

```
cat summary.txt  | awk '{print $1}' > u
for s in `cat u`; do var1=`find . -name *${s%.[0-9]}*`; echo $s$'\t'$var1 >> ncbiID_to_file; done
```

* it turns out that this is a fairly comprehensive list 25 genomes were dropped because they were not found in the MetaData

```
grep --perl-regexp "\t$" ncbiID_to_file | wc
      25      25     324
grep --perl-regexp "\t$" ncbiID_to_file | wc
NC_003911.11  
AC_000091.1	
NC_009353.1	
NC_009444.1	
NC_010332.1	
NS_000190.1	
NC_011980.1	
NS_000196.1	
NS_000197.2	
NC_012627.1	
NC_012629.1	
NC_012630.1	
NC_012915.1	
NC_013416.1	
NC_013438.1	
NC_013597.1	
NC_013784.1	
NC_013785.1	
NC_013786.1	
NC_013787.1	
NC_013788.1	
NC_014629.1	
NC_015557.1	
NC_015587.1
```

* this leaves us with 2617 genomes with annotation of 2642, storing these in `ncbiID_to_file.txt`

```
grep --perl-regexp ".*fna" ncbiID_to_file2 | wc
    2617    5234  197129
wc ncbiID_to_file2
    2642    5259  197453 ncbiID_to_file2
grep --perl-regexp ".*fna" ncbiID_to_file2 > ncbiID_to_file.txt
```

## Simulating Genomes

Here, I create two simputated metagenomes consisting of 10,000bp contigs randomly generated from 2713 NCBI genomes (Downloaded March 15 2014 from <ftp://ftp.ncbi.nlm.nih.gov/pub/taxonomy/>) using the python script `subsample_ncbi.py`. The first simulation is a smaller sample of 100 genomes, sampling 10 random reads from each genome. The second simulation is larger, sampling a random subset of 2000 genomes, sampling 10 random reads from each. Each simulation had its ORFs predicted and annotated against the RefSeq database using the MetaPathways pipeline. The original source of each contig was predicted from the taxonomic annotations ascribed to each LCA-predicted ORF using three different methods: $LCA^2$, $Simple Majority$, and our information-theoretic $LCA*$. $LCA^2$ simply applies the LCA algorithm again to the set of contig taxonomies. The $Simple Majority$ method ascribes the taxonomy of the contig to the taxonomy that has the greatest majority. Our $LCA*$ method applies the information-theoetic result and algorithm previously described with a majority theshold set to the default majority ($\alpha=0.5$).

We evaluated the performance of these predictions using two taxonomic distances on the NCBI taxonomy tree. First a simple walk on the NCBI Taxonomy Hierarchy from the observed predicted taxonomy to the original expected taxonomy. The second is a weighted taxonomic distance that weightes each edge proportional to $\frac{1}{d}$ where $d$ is the depth of the edge in the tree. See Hanson *et al.* (2014) for more details. The NCBI Taxonomy Hierarchy was modified with the additional *prokaryotes* node as a parent of *Bacteria* and *Archaea* nodes.

### Creating Simulated Metagenomes

Wrote my own script `subsample_ncbi.py` quickly sample from a collection of fasta files specified in our `ncbiID_to_file.txt` that we created above.

The following example creates sub-sequences of length 10,000, sampling 10 sequences per file on a random subset of 100 (the random number generator is seeded so that results can be reproducible)

```
python subsample_ncbi.py -i ncbiID_to_file.txt -l 10000 -s 10 -n 100 -o lca_star_test1.fna
```

* lca_star_test1.fna is going to be our first test, we ran it through MetaPathways with the standard settings:

```
Run Date : 2014-03-15 
Nucleotide Quality Control parameters
  min length  180
ORF prediction parameters
  min length	60
  algorithm	prodigal
Amino acid quality control and annotation parameters
  min bit score	20
  min seq length	60
  annotation reference dbs	RefSeq_complete_nr_v62_dec2013
  min BSR	0.4
  max evalue	0.000001
Pathway Tools parameters
  taxonomic pruning 	no
rRNA search/match parameters
  min identity	20
  max evalue	0.000001
  rRNA reference dbs	GREENGENES_gg16S-2013-09-12
```

### C

* `create_analysis_table.py` summarizes results into a Long Table format that can be processed by R:

```
python create_analysis_table.py -i lca_star_test/metapathways_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt --ft_table lca_star_test/metapathways_out/lca_star_test1/results/annotation_table/functional_and_taxonomic_table.txt --sum_table lca_star_test/ncbi_bacteria_genomes/ncbi/summary.txt --ncbi_tree ~/Dropbox/utilities/tree_update/data/pipeline/ncbi_taxonomy_tree.txt -o test.out.txt
```

* example output:

```
shebop:lca_star_test nielsh$ head first_test.txt 
contig	real	taxa	method	dist
lca_star_test1_378	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	Desulfomonas Moore et al. 1976	LCA_Star	3
lca_star_test1_378	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	all	Majority	11
lca_star_test1_378	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	biota	LCA_Squared	10
lca_star_test1_379	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	Vibrio cholinicus	LCA_Star	2
lca_star_test1_379	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	Vibrio cholinicus	Majority	2
lca_star_test1_379	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	Desulfomonas Moore et al. 1976	LCA_Squared	3
lca_star_test1_372	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	Desulfomonas Moore et al. 1976	LCA_Star	3
lca_star_test1_372	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	all	Majority	11
lca_star_test1_372	Desulfovibrio desulfuricans subsp. desulfuricans str. ATCC 27774	'Desulfovibrionaceae'	LCA_Squared	4
```

* a function `lcastar.get_distance(taxon, real)` was added to the LCAStar.py class to calculate a simple distance between taxa.
* TODO: come up with a distance of 'divergence' where we penalize predicted lineages that are outside of the true lineage

### Preliminary statistics

First, load required libraries. 


```r
library(ggplot2)
```

#### First test

Here compared our first implementation (March 16 2014) of LCA* using the `lca_star_test1.fasta` above (10,000, sampling 10 sequences per file on a random subset of 100):













