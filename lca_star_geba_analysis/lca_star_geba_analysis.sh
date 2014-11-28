#!/bin/bash

projects_root=/Users/nielshanson/Dropbox/projects

# geba project
input=${projects_root}/LCAStar_data/assemblies/ncbi_names/sag_ids_to_ncbi.txt
mp_out=${projects_root}/LCAStar_data/assemblies/GEBA_SAGs_out/
output=${projects_root}/LCAStar/lca_star_geba_analysis/
ncbi_tree_dir=/Users/nielshanson/Desktop/MetaPathways_DBs/ncbi_tree
while read sag taxa
do
    blast_result=$mp_out/${sag}/blast_results/${sag}.RefSeq_nr_v62_MDM_SAGs_dec2013.LASTout.parsed.txt
    mapping_file=$mp_out/${sag}/preprocessed/${sag}.mapping.txt
    echo python2.7 ${projects_root}/LCAStar/Compute_LCAStar.py \
           -i $blast_result \
           -m $mapping_file \
           --ncbi_tree ${ncbi_tree_dir}/ncbi_taxonomy_tree.txt \
           --ncbi_megan_map ${ncbi_tree_dir}/ncbi.map \
           --orf_summary besthit \
           --sample_taxa_ref \"${taxa}\" \
           --all_methods \
		   -o $output/${sag}_lcastar.txt
done < $input

# test1 and test2
# creating contig_maps
# python contig_to_taxa_map.py -i metapathways_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt -n ncbi_bacteria_genomes/ncbi/summary.txt -o lca_star_test1_contigmap.txt
# python contig_to_taxa_map.py -i metapathways_out/lca_star_test1/preprocessed/lca_star_test2.mapping.txt -n ncbi_bacteria_genomes/ncbi/summary.txt -o lca_star_test2_contigmap.txt
mp_out=/Users/nielshanson/Dropbox/projects/LCAStar_data/lca_star_test/metapathways_out

for sample in test1 test2
do
    blast_result=$mp_out/lca_star_${sample}/blast_results/lca_star_${sample}.refseq.lastout.parsed.txt
    mapping_file=$mp_out/lca_star_${sample}/preprocessed/lca_star_${sample}.mapping.txt
    contig_ref=/Users/nielshanson/Dropbox/projects/LCAStar_data/lca_star_test/lca_star_${sample}_contigmap.txt
    output=${projects_root}/LCAStar/lca_star_geba_analysis/
    echo python2.7 ${projects_root}/LCAStar/Compute_LCAStar.py \
               -i $blast_result \
               -m $mapping_file \
               --ncbi_tree ${ncbi_tree_dir}/ncbi_taxonomy_tree.txt \
               --ncbi_megan_map ${ncbi_tree_dir}/ncbi.map \
               --orf_summary lca \
               --contig_taxa_ref $contig_ref \
               --all_methods \
               -o $output/${sample}_lcastar.txt
done
