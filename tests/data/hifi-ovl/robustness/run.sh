#! /bin/bash

# Extracts the reads from the dataset.

seqdb=/pbi/dept/secondary/siv/testdata/ipa/ipa2/dros_test_data/haplotype_tagged_38x_19kb_hifi_reads/reads.seqdb

for f in reads.robustness.*.txt
do
    echo $f
    pancake seqfetch ${f%.*}.fasta ${f} ${seqdb}
done
