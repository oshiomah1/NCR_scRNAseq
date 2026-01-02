#!/bin/bash

# Define the output file
output_file="/share/hennlab/projects/NCR_scRNAseq/results/demultiplex/freemuxlet/1000g_remove.csv"

# Add header from the first file, then append all files, skipping headers after the first file
{
  head -n 1 "/share/hennlab/projects/CAAPA2_functional_annotation/allele_frequencies_and_scores/data/pop_filtered/pop_list/carribean/Barbados1000G.csv"
  for file in /share/hennlab/projects/CAAPA2_functional_annotation/allele_frequencies_and_scores/data/pop_filtered/pop_list/carribean/Barbados1000G.csv \
              /share/hennlab/projects/CAAPA2_functional_annotation/allele_frequencies_and_scores/data/pop_filtered/pop_list/carribean/PuertoRican1000G.csv \
              /share/hennlab/projects/CAAPA2_functional_annotation/allele_frequencies_and_scores/data/pop_filtered/pop_list/WestAf/Yoruba1000G.csv
  do
    tail -n +2 "$file"
  done
} > "$output_file"
