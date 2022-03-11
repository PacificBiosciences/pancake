These were generated like so:

```
${BIN_DIR}/pancake seeddb -k 30 -w 80 ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1a

${BIN_DIR}/pancake seeddb -k 30 -w 80 --split-blocks ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1b
```

The `test-2a.seeddb` and `test-2b.seeddb` are the same as `test-1a.seeddb`, but the number of bytes for each SeedsLine was manually changed to an invalid value for `test-2a`, while for `test-2b` the number of seeds was changed to an invalid value.
