# ecoli fasta from https://www.ncbi.nlm.nih.gov/nuccore/U00096.2?report=fasta&log$=seqview&format=text

CPU=14

# MP3_REF=/home/tony/workspace/resources/mp3db
# metapathways run \
#     -t $CPU \
#     -i ./cache/ecoli.fna \
#     -o ./cache/ecoli_mp3 \
#     --COMPUTE_TPM skip \
#     -d $MP3_REF

NR_DMND=/home/tony/workspace/resources/nr/nr.dmnd
diamond blastp \
    --threads $CPU \
    --outfmt 6 qseqid sseqid stitle pident bitscore \
    --db $NR_DMND \
    --query ./cache/ecoli_mp3/ecoli/orf_prediction/ecoli.faa \
    --out ./cache/ecoli_nr.tsv
