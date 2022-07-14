#!/bin/bash
set -ue

region_fasta_file=$1
background_region_fasta_file=$2
motif_file=$3
output_folder=$4

sea $region_fasta_file $motif_file -o $output_folder -n $background_region_fasta_file --qvalue