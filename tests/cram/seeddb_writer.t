Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ ${BIN_DIR}/pancake seeddb -k 20 -w 80 ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1a
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1a.seeddb test-1a.seeddb
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1a.seeddb.0.seeds test-1a.seeddb.0.seeds
  > ls -1 test-1a.seeddb*
  test-1a.seeddb
  test-1a.seeddb.0.seeds

Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ ${BIN_DIR}/pancake seeddb -k 20 -w 80 --split-blocks ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1b
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1b.seeddb test-1b.seeddb
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1b.seeddb.0.seeds test-1b.seeddb.0.seeds
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1b.seeddb.1.seeds test-1b.seeddb.1.seeds
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1b.seeddb.2.seeds test-1b.seeddb.2.seeds
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1b.seeddb.3.seeds test-1b.seeddb.3.seeds
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-1b.seeddb.4.seeds test-1b.seeddb.4.seeds
  > ls -1 test-1b.seeddb*
  test-1b.seeddb
  test-1b.seeddb.0.seeds
  test-1b.seeddb.1.seeds
  test-1b.seeddb.2.seeds
  test-1b.seeddb.3.seeds
  test-1b.seeddb.4.seeds

Build the SeedDB from a SeqDB which is split into multiple blocks, and split the seeds into multiple files.
  $ rm -f test-11*
  > ${BIN_DIR}/pancake seeddb --split-blocks -k 20 -w 80 ${TEST_DATA_DIR}/seqdb-writer/test-11-split-blocks.seqdb test-11b-split-blocks
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-11b-split-blocks.seeddb test-11b-split-blocks.seeddb
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-11b-split-blocks.seeddb.0.seeds test-11b-split-blocks.seeddb.0.seeds
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-11b-split-blocks.seeddb.1.seeds test-11b-split-blocks.seeddb.1.seeds
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-11b-split-blocks.seeddb.2.seeds test-11b-split-blocks.seeddb.2.seeds
  > ls -1 test-11b-split-blocks.seeddb*
  test-11b-split-blocks.seeddb
  test-11b-split-blocks.seeddb.0.seeds
  test-11b-split-blocks.seeddb.1.seeds
  test-11b-split-blocks.seeddb.2.seeds

Build the SeedDB from a SeqDB which is split into multiple blocks, but keeps the seeds in a single file.
  $ rm -f test-11*
  > ${BIN_DIR}/pancake seeddb -k 20 -w 80 ${TEST_DATA_DIR}/seqdb-writer/test-11-split-blocks.seqdb test-11a-single-file
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-11a-single-file.seeddb test-11a-single-file.seeddb
  > diff ${TEST_DATA_DIR}/seeddb-writer/test-11a-single-file.seeddb.0.seeds test-11a-single-file.seeddb.0.seeds
  > ls -1 test-11a-single-file.seeddb*
  test-11a-single-file.seeddb
  test-11a-single-file.seeddb.0.seeds
