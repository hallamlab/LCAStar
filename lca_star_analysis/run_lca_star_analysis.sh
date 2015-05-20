#!/bin/bash

# lca_star_geba_analysis.sh
# Simple script to run LCAStar on the Small and Large simulated metagenome contigs
# as well as the GEBA Single-cell Amplified Genomes (SAGs).

# BASE PATH - Update to re-run
lca_star=/Users/nielshanson/Dropbox/projects/LCAStar 
lca_star_data=${lca_star}/lca_star_data

mp_out=${lca_star_data}/assemblies/GEBA_SAGs_out/



# Creating contig_maps for the Small (test1) and Large (test2) simulations

lca_star_test=${lca_star_data}/mp_data/lca_star_simulations
#echo python2.7 contig_to_taxa_map.py -i ${lca_star_test}/mp_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt \
#                                -n ${lca_star_data}/ncbi_bacteria_genomes/ncbi/summary.txt \
#                                -o ${lca_star_test}/lca_star_test1_contigmap.txt
#echo python2.7 contig_to_taxa_map.py -i ${lca_star_test}/mp_out/lca_star_test2/preprocessed/lca_star_test2.mapping.txt \
#                                -n ${lca_star_data}/ncbi_bacteria_genomes/ncbi/summary.txt \
#                                -o ${lca_star_test}/lca_star_test2_contigmap.txt

mp_out=${lca_star_test}/mp_out
ncbi_tree_dir=${lca_star}/resources
lca_star_out=${lca_star}/lca_star_analysis/lca_star_results/

# Run Small and Large simluations
for sample in test1 test2
do
    blast_result=${mp_out}/lca_star_${sample}/blast_results/lca_star_${sample}.refseq.lastout.parsed.txt
    mapping_file=${mp_out}/lca_star_${sample}/preprocessed/lca_star_${sample}.mapping.txt
    contig_ref=${lca_star_test}/lca_star_${sample}_contigmap.txt
    output=${lca_star}/lca_star_geba_analysis/
    echo python2.7 ${lca_star}/Compute_LCAStar.py \
               -i $blast_result \
               -m $mapping_file \
               --ncbi_tree ${ncbi_tree_dir}/ncbi_taxonomy_tree.txt \
               --ncbi_megan_map ${ncbi_tree_dir}/ncbi.map \
               --orf_summary lca \
               --contig_taxa_ref ${contig_ref} \
               --all_methods \
               -o ${lca_star_out}/${sample}_lcastar.txt
done


# Run GEBA MDM SAG analysis

# input=${lca_star_data}/assemblies/ncbi_names/sag_ids_to_ncbi.txt
# output=${lca_star}/lca_star_geba_analysis/

# while read sag taxa
# do
#     blast_result=$mp_out/${sag}/blast_results/${sag}.RefSeq_nr_v62_MDM_SAGs_dec2013.LASTout.parsed.txt
#     mapping_file=$mp_out/${sag}/preprocessed/${sag}.mapping.txt
#     echo python2.7 ${projects_root}/LCAStar/Compute_LCAStar.py \
#            -i $blast_result \
#            -m $mapping_file \
#            --ncbi_tree ${ncbi_tree_dir}/ncbi_taxonomy_tree.txt \
#            --ncbi_megan_map ${ncbi_tree_dir}/ncbi.map \
#            --orf_summary besthit \
#            --sample_taxa_ref \"${taxa}\" \
#            --all_methods \
#        -o $output/${sag}_lcastar.txt
# done < $input

