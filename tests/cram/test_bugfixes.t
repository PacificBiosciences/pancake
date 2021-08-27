############################
Labeling of overlaps. There was a bug where one of the overlaps used to be marked as 3 prime instead of 5 prime.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/bugfixes/overlap_labelling/reads.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --skip-sym --write-rev --num-threads 1 reads reads 0 0 0
  pat_m64011_190605_003147/160434972/ccs pat_m64011_190605_003147/32900288/ccs -11692 99.78 0 0 11718 13995 1 0 11737 13966 5
  pat_m64011_190605_003147/32900288/ccs pat_m64011_190605_003147/160434972/ccs -11692 99.78 0 0 11737 13966 1 0 11718 13995 5

Labeling of overlaps, same as before but write IDs instead of headers.
  $ ${BIN_DIR}/pancake seqdb reads ${PROJECT_DIR}/test-data/bugfixes/overlap_labelling/reads.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --skip-sym --write-rev --write-ids --num-threads 1 reads reads 0 0 0
  000000001 000000000 -11692 99.78 0 0 11718 13995 1 0 11737 13966 5
  000000000 000000001 -11692 99.78 0 0 11737 13966 1 0 11718 13995 5

############################
### Empty SeqDB as input ###
############################
SeqDB test for empty input.
  $ rm -rf out; mkdir -p out
  > touch out/in.fasta
  > ${BIN_DIR}/pancake seqdb out/out out/in.fasta

SeedDB test for empty SeqDB.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seeddb ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb out/out

Ovl-hifi test for empty SeqDB.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake ovl-hifi ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs 0 0 0

DBFilter test for empty SeqDB.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake dbfilter ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs out/out

Seqfetch test for empty SeqDB.
  $ rm -rf out; mkdir -p out
  > touch out/in.txt
  > ${BIN_DIR}/pancake seqfetch out/out.fasta out/in.txt ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb

Seqdb-dump test for empty SeqDB.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb out/out.fasta

Seqdb-info test for empty SeqDB.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb --human | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	0.00	0	0.00	0.00	0.00	0.00	0.00	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0

Composite test for multiple tools, from an empty input FASTA.
  $ rm -rf out; mkdir -p out
  > touch out/in.fasta
  > ${BIN_DIR}/pancake seqdb out/out out/in.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0
  > ${BIN_DIR}/pancake dbfilter out/out out/out.filtered
  > touch out/to_fetch.txt
  > ${BIN_DIR}/pancake seqfetch out/out.fetched.fasta out/to_fetch.txt out/out.seqdb
  > wc -l out/out.fetched.fasta | awk '{ print $1 }'
  > ${BIN_DIR}/pancake seqdb-dump out/out.seqdb out/out.all.fasta
  > diff out/in.fasta out/out.all.fasta
  > ${BIN_DIR}/pancake seqdb-info out/out.seqdb --human | cut -f 2-
  0
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	0.00	0	0.00	0.00	0.00	0.00	0.00	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0
