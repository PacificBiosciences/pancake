The dataset can be recreated with:
```
pancake seqdb --buffer-size 1024 --block-size 0 out in.fasta
...
pancake seqdb --compression 0 --block-size 0.03 test-7-uncompressed-2blocks in.fasta
pancake seqdb --compression 1 --block-size 0.03 test-8-compressed-2blocks in.fasta
```
