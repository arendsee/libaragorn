#!/usr/bin/env bash

./test data/sample-trna_sequence_cmp_arc_1.fasta tRNA > test-trna_sequence_cmp_arc_1.txt
echo $(grep -c '>' data/sample-trna_sequence_cmp_arc_1.fasta) $(grep -c "seq = " test-trna_sequence_cmp_arc_1.txt)

./test data/sample-trna_sequence_cmp_bac_1.fasta tRNA > test-trna_sequence_cmp_bac_1.txt
echo $(grep -c '>' data/sample-trna_sequence_cmp_bac_1.fasta) $(grep -c "seq = " test-trna_sequence_cmp_bac_1.txt)

./test data/sample-trna_sequence_fungi_1.fasta tRNA > test-trna_sequence_fungi_1.txt
echo $(grep -c '>' data/sample-trna_sequence_fungi_1.fasta) $(grep -c "seq = " test-trna_sequence_fungi_1.txt)

./test data/sample-trna_sequence_phage_1.fasta tRNA > test-trna_sequence_phage_1.txt
echo $(grep -c '>' data/sample-trna_sequence_phage_1.fasta) $(grep -c "seq = " test-trna_sequence_phage_1.txt)

./test data/sample-trna_sequence_plant_1.fasta tRNA > test-trna_sequence_plant_1.txt
echo $(grep -c '>' data/sample-trna_sequence_plant_1.fasta) $(grep -c "seq = " test-trna_sequence_plant_1.txt)

./test data/sample-trna_sequence_plasmid_1.fasta tRNA > test-trna_sequence_plasmid_1.txt
echo $(grep -c '>' data/sample-trna_sequence_plasmid_1.fasta) $(grep -c "seq = " test-trna_sequence_plasmid_1.txt)

./test data/sample-trna_sequence_virus_1.fasta tRNA > test-trna_sequence_virus_1.txt
echo $(grep -c '>' data/sample-trna_sequence_virus_1.fasta) $(grep -c "seq = " test-trna_sequence_virus_1.txt)

./test data/Mycoplasmoides_pneumoniae.fasta tRNA > test-Mycoplasmoides_pneumoniae.txt
echo $(grep -c "seq = " test-Mycoplasmoides_pneumoniae.txt)

./test data/sample-bac-intronic.fa tRNA > test-bac-intronic.txt
echo $(grep -c '>' data/sample-bac-intronic.fa) $(grep -c "seq = " test-bac-intronic.txt)

./test data/sample-euk-intronic.fa tRNA > test-euk-intronic.txt
echo $(grep -c '>' data/sample-euk-intronic.fa) $(grep -c "seq = " test-euk-intronic.txt)


echo "Testing mtRNA"

./test data/rat-mito.fna mtRNA > z


# If a file of reference results is given, then compare to current results
if [[ -d gold ]]
then
  for j in test-*; do diff $j gold/$j; done
fi
