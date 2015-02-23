#!/bin/bash

# lca_star_geba_analysis.sh
# Simple script to run LCAStar on the Small and Large simulated metagenome contigs
# as well as the GEBA Single-cell Amplified Genomes (SAGs).

projects_root=/Users/nielshanson/Dropbox/projects
input=${projects_root}/LCAStar_data/assemblies/ncbi_names/sag_ids_to_ncbi.txt
mp_out=${projects_root}/LCAStar_data/assemblies/GEBA_SAGs_out/
output=${projects_root}/LCAStar/lca_star_geba_analysis/
ncbi_tree_dir=${projects_root}/LCAStar/resources
mp_out=${lca_star_test}/metapathways_out

# test1 and test2
# creating contig_maps

lca_star_test=${projects_root}/LCAStar_data/lca_star_test
python2.7 contig_to_taxa_map.py -i ${lca_star_test}/metapathways_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt \
                             -n ${lca_star_test}/ncbi_bacteria_genomes/ncbi/summary.txt \
                             -o ${lca_star_test}/resources/lca_star_test1_contigmap.txt
python2.7 contig_to_taxa_map.py -i ${lca_star_test}/metapathways_out/lca_star_test2/preprocessed/lca_star_test2.mapping.txt \
                                -n ${lca_star_test}/ncbi_bacteria_genomes/ncbi/summary.txt \
                                -o ${lca_star_test}/resources/lca_star_test2_contigmap.txt

for sample in test1 test2
do
    blast_result=${mp_out}/lca_star_${sample}/blast_results/lca_star_${sample}.refseq.lastout.parsed.txt
    mapping_file=${mp_out}/lca_star_${sample}/preprocessed/lca_star_${sample}.mapping.txt
    contig_ref=${projects_root}/LCAStar_data/lca_star_test/lca_star_${sample}_contigmap.txt
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

