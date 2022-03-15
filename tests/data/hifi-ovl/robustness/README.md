The data in this folder attempts to cover the widest spread of potential overlap types that can happen:
- Regular plain (linear) reads.
- Tandem repeats.
- Low-complexity reads.
- Reads with various artefacts.
- Reads with suprious hits when mapping. Without proper refining of seed hits, this can cause the overlap to have a bad alignment and be filtered out. (Example: `reads.robustness.spurious_hits.1.*.fasta`.)

All the data in this folder is extracted from the Drosophila dataset.
The dataset can be found here:
/pbi/dept/secondary/siv/testdata/ipa/ipa2/dros_test_data/haplotype_tagged_38x_19kb_hifi_reads/reads.seqdb

The `run.sh` script is intended only for documentation/reproducibility - it demonstrates how the reads were extracted.
