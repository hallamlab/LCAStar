LCAStar: an entropy-based measure for taxonomic assignment within assembled metagenomes
=======

Niels W. Hanson, Kishori M. Konwar, Steven J. Hallam

![lca_star_logo.png](lca_star_logo.png)

## Abstract

A perennial problem in the analyses of large meta'omic datasets is the taxonomic classification of unknown reads or assembled contigs to their likely taxa of origin. Although the assembly of metagenomic samples has its difficulties, once contigs are found it is often important to classify them to a taxonomy based on their ORF annotations. The popular Lowest Common Ancestor (LCA) algorithm addresses a similar problem with ORF annotations, and it is intuitive to apply the same taxonomic annotation procedure to the annotation of contigs, a procedure we call LCA2. Inspired by Information and Voting Theory we developed an alternative statistics LCA\* by viewing the taxonomic classification problem as an election among the different taxonomic annotations, and gen- eralize an algorithm to obtain a sufficiently strong majority α-majority while respecting the entropy of the taxonomic distribution and phylogeny tree-structure of the NCBI Taxonomic Database. Further, using results from order and supremacy statistics, we formulate a likelihood-ratio hypothesis test and p-value for testing the supremacy of the final reported taxonomy. In simulated metage- nomic config experiments, we emperically demonstrate that voting-based methods, majority vote and LCA\*, are significantly more accurate than LCA2, and that in many cases LCA\* is superior to the simple majority vote procedure. LCA\* and its statistical tests have been implemented as a stand-alone Python library, and have been integrated into the latest release of the [MetaPathways pipeline](https://github.com/hallamlab/metapathways2).

## Installation

LCA\* is released as as Python library, requiring Python 2.6 or greater. More installation and useage information can be found on the wiki.

## Contents

* [Compute_LCAStar.py](Compute_LCAStar.py): Driver script for running LCAStar.py

    * Usage: 
    
    ```
    python Compute_LCAStar.py -i blast_results/refseq.*.parsed.txt \
                              -m preprocessed/*.mapping.txt \
                              --ncbi_tree resources/ncbi_taxonomy_tree.txt \
                              --ncbi_megan_map resources/ncbi.map \
                              -a \
                              -v \
                              --contig_taxa_ref ...contigmap.txt \
                              -o LCAStar.output.txt
    ```
    
    where,
    
        * `-i`: is a MetaPathways `parsed.txt` annotation file
        * `-m`: is a MetaPathways mapping file `.mapping.txt`
        * `--ncbi_tree`: the MetaPathways `ncbi_taxonomy_tree.txt`
        * `-a`: computes all methods Majority, LCAStar, and LCA^2
        * `-v`: verbose mode
        * `--contig_taxa_ref`: calculat
        * `-o`: output text file


