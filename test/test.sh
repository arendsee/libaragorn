#!/usr/bin/env bash

make

./test data/sample-trna_sequence_cmp_arc_1.fasta > test-trna_sequence_cmp_arc_1.txt
echo $(grep -c '>' data/sample-trna_sequence_cmp_arc_1.fasta) $(wc -l test-trna_sequence_cmp_arc_1.txt)

./test data/sample-trna_sequence_cmp_bac_1.fasta > test-trna_sequence_cmp_bac_1.txt
echo $(grep -c '>' data/sample-trna_sequence_cmp_bac_1.fasta) $(wc -l test-trna_sequence_cmp_bac_1.txt)

./test data/sample-trna_sequence_fungi_1.fasta > test-trna_sequence_fungi_1.txt
echo $(grep -c '>' data/sample-trna_sequence_fungi_1.fasta) $(wc -l test-trna_sequence_fungi_1.txt)

./test data/sample-trna_sequence_phage_1.fasta > test-trna_sequence_phage_1.txt
echo $(grep -c '>' data/sample-trna_sequence_phage_1.fasta) $(wc -l test-trna_sequence_phage_1.txt)

./test data/sample-trna_sequence_plant_1.fasta > test-trna_sequence_plant_1.txt
echo $(grep -c '>' data/sample-trna_sequence_plant_1.fasta) $(wc -l test-trna_sequence_plant_1.txt)

./test data/sample-trna_sequence_plasmid_1.fasta > test-trna_sequence_plasmid_1.txt
echo $(grep -c '>' data/sample-trna_sequence_plasmid_1.fasta) $(wc -l test-trna_sequence_plasmid_1.txt)

./test data/sample-trna_sequence_virus_1.fasta > test-trna_sequence_virus_1.txt
echo $(grep -c '>' data/sample-trna_sequence_virus_1.fasta) $(wc -l test-trna_sequence_virus_1.txt)

./test data/Mycoplasmoides_pneumoniae.fasta > test-Mycoplasmoides_pneumoniae.txt
echo $(wc -l test-Mycoplasmoides_pneumoniae.txt)
