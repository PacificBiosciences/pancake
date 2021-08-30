Tests dumping of the entire DB to FASTA.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta out.fasta

Tests dumping of the entire DB to FASTA. Write sequence IDs instead of names.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta --write-ids
  > awk 'BEGIN {n = 0} { if (substr($1,1,1) == ">") { print ">00000000"n; ++n } else { print } }' ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta > expected.fasta
  > diff expected.fasta out.fasta

Tests dumping of the entire DB to stdout.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb - > out.stdout.fasta
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta out.stdout.fasta

Test dumping of a particular block ID.
Note - the Awk command here concatenates the multiline FASTA into a single line because Samtools automatically wraps.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta --block-id 2
  > samtools faidx ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" | awk 'BEGIN { n = 0} { c = substr($1, 1, 1); if (c == ">") { if (n > 0) { printf("\n") }; printf("%s\n", $0)} else { printf("%s", $0) } n = 1} END { printf("\n") }' > expected.fasta
  > diff expected.fasta out.fasta

Test dumping of the first block only.
Note - the Awk command here concatenates the multiline FASTA into a single line because Samtools automatically wraps.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta --block-id 0
  > samtools faidx ${PROJECT_DIR}/test-data/seqdb-writer/in.fasta "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852" | awk 'BEGIN { n = 0} { c = substr($1, 1, 1); if (c == ">") { if (n > 0) { printf("\n") }; printf("%s\n", $0)} else { printf("%s", $0) } n = 1} END { printf("\n") }' > expected.fasta
  > diff expected.fasta out.fasta

Try dumping a block out of bounds.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb out.fasta --block-id 10 2>&1 | sed 's/.*pancake //g'
  seqdb-dump ERROR: Specified block ID is too large, numBlocks = 5.

Test fetching the homopolymer compressed sequences.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/seqdb-writer/test-10a-compressed-in-2-small.seqdb - --use-hpc > out.fasta
  > diff ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.hpc_compressed.fasta out.fasta
