# Run some regression tests

python Compute_LCAStar.py -i /Users/nielshanson/Dropbox/projects/LCAStar_data/lca_star_test/metapathways_out/lca_star_test1/blast_results/lca_star_test1.refseq.lastout.parsed.txt -m /Users/nielshanson/Dropbox/projects/LCAStar_data/lca_star_test/metapathways_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt --ncbi_tree /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi_taxonomy_tree.txt --ncbi_megan_map /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi.map -a -v --contig_taxa_ref /Users/nielshanson/Dropbox/projects/LCAStar_data/lca_star_test/resources/lca_star_test1_contigmap.txt 


# Run new simulations
python Compute_LCAStar.py -i /Users/nielshanson/metapathways/metapathways2/output/lca_star_test1/blast_results//lca_star_test1.refseq_tests_removed_03-04-2016.LASTout.parsed.txt \
                          -m /Users/nielshanson/Dropbox/projects/LCAStar/lca_star_data/mp_data/lca_star_simulations/mp_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt \
                          --ncbi_tree /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi_taxonomy_tree.txt \
                          --ncbi_megan_map /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi.map \
                          --orf_summary lca \
                          --all_methods \
                          -a \
                          --contig_taxa_ref /Users/nielshanson/Dropbox/projects/LCAStar/lca_star_data/mp_data/lca_star_simulations/lca_star_test1_contigmap.txt \
                          -o /Users/nielshanson/Desktop/lca_star_test1_removed_best_03-04-2016_lcastar.txt

python Compute_LCAStar.py -i /Users/nielshanson/metapathways/metapathways2/output/lca_star_test2/blast_results//lca_star_test2.refseq_tests_removed_03-04-2016.LASTout.parsed.txt \
                          -m /Users/nielshanson/Dropbox/projects/LCAStar/lca_star_data/mp_data/lca_star_simulations/mp_out/lca_star_test2/preprocessed/lca_star_test2.mapping.txt \
                          --ncbi_tree /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi_taxonomy_tree.txt \
                          --ncbi_megan_map /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi.map \
                          --orf_summary lca \
                          --all_methods \
                          -a \
                          --contig_taxa_ref /Users/nielshanson/Dropbox/projects/LCAStar/lca_star_data/mp_data/lca_star_simulations/lca_star_test2_contigmap.txt \
                          -o /Users/nielshanson/Desktop/lca_star_test2_removed_best_03-04-2016_lcastar.txt

# python ${lca_star}/Compute_LCAStar.py \
#          -i $blast_result \
#          -m $mapping_file \
#          --ncbi_tree ${ncbi_tree_dir}/ncbi_taxonomy_tree.txt \
#          --ncbi_megan_map ${ncbi_tree_dir}/ncbi.map \
#          --orf_summary lca \
#          --contig_taxa_ref ${contig_ref} \
#          --all_methods \
#          -o ${lca_star_out}/${sample}_lcastar.txt

python Compute_LCAStar.py -i /Users/nielshanson/metapathways/metapathways2/output/lca_star_test1/blast_results//lca_star_test1.refseq_tests_removed_03-04-2016.LASTout.parsed.txt                           -m /Users/nielshanson/Dropbox/projects/LCAStar/lca_star_data/mp_data/lca_star_simulations/mp_out/lca_star_test1/preprocessed/lca_star_test1.mapping.txt                           --ncbi_tree /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi_taxonomy_tree.txt                           --ncbi_megan_map /Users/nielshanson/Dropbox/projects/LCAStar/resources/ncbi.map                           --orf_summary                           --all_methods                           -a                           --contig_taxa_ref /Users/nielshanson/Dropbox/projects/LCAStar/lca_star_data/mp_data/lca_star_simulations/lca_star_test1_contigmap.txt                           -o ~/Desktop/lca_star_test1_removed_03-04-2016_lcastar.txt