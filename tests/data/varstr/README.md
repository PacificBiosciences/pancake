The pair of reads is taken from the synthetic dataset ivan-200k-t1.

```
function pair1 {
    local out_prefix=reads1.weird_masking_issue
    echo "000000551" > ${out_prefix}.txt
    echo "000000571" >> ${out_prefix}.txt
    SEQ_DB_PREFIX="../RUN/02-build_db/reads"
    pancake seqfetch --write-ids --alias ${SEQ_DB_PREFIX}.seqdb ${out_prefix}.fasta ${out_prefix}.txt ${SEQ_DB_PREFIX}.seqdb
}
```

The reads3.masking_and_identity.fasta file consists of a couple of reads from the Murder Hornet dataset. On this dataset, I noticed that the resulting overlap did not have 100% identity.
