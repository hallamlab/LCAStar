#!/bin/bash

input=/Users/nielshanson/Dropbox/projects/LCAStar_data/assemblies/ncbi_names/sag_ids_to_ncbi.txt
mp_out=/Users/nielshanson/Dropbox/projects/LCAStar_data/assemblies/GEBA_SAGs_out/
output=/Users/nielshanson/Dropbox/projects/LCAStar/lca_star_geba_analysis/
while read sag taxa
do
    blast_result=$mp_out/$sag/blast_results/$sag.RefSeq_nr_v62_MDM_SAGs_dec2013.LASTout.parsed.txt
    mapping_file=$mp_out/$sag/preprocessed/$sag.mapping.txt
    echo python /Users/nielshanson/Dropbox/projects/LCAStar/Compute_LCAStar.py \
           -i $blast_result \
           -m $mapping_file \
           --ncbi_tree /Users/nielshanson/Downloads/MetaPathways_DBs/ncbi_tree/ncbi_taxonomy_tree.txt \
           --ncbi_megan_map /Users/nielshanson/Downloads/MetaPathways_DBs/ncbi_tree/ncbi.map \
           --orf_summary besthit \
           --sample_taxa_ref \"${taxa}\" \
           --all_methods \
		   -o $output/${sag}_lcastar.txt
done < $input


