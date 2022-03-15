##############################
### Filtering with a list. ###
##############################
Wrong filter, should throw.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --filter-type something_wrong ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered 2>&1 | sed 's/^.*pancake //g'
  dbfilter ERROR: Unknown filter type: 'something_wrong'.

No filtering applied.
Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  test-1.filtered.seqdb

No filtering, even though a filter list is specified, but the filter type is none.
Filter out two sequences using a blacklist. The actual data files will be reused from the input DB.
Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ rm -f test-1.filtered*
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" > filterlist.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> filterlist.txt
  $ ${BIN_DIR}/pancake dbfilter --block-size 0 --filter-list filterlist.txt --filter-type none ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  test-1.filtered.seqdb

Filter out two sequences using a blacklist. The actual data files will be reused from the input DB.
Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ rm -f test-1.filtered*
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" > filterlist.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> filterlist.txt
  $ ${BIN_DIR}/pancake dbfilter --block-size 0 --filter-list filterlist.txt --filter-type blacklist ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  V	0.1.0
  C	1
  F	0	test-1.seqdb.0.seq	1	1463	5852
  F	1	test-1.seqdb.1.seq	1	2996	11983
  F	2	test-1.seqdb.2.seq	1	6073	24292
  F	3	test-1.seqdb.3.seq	1	1277	5105
  F	4	test-1.seqdb.4.seq	1	4751	19001
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	2	0	6073	24292	1	0	24292
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	3	0	1277	5105	1	0	5105
  B	0	0	1	1463	5852
  B	1	1	2	6073	24292
  B	2	2	3	1277	5105
  test-1.filtered.seqdb

Filter out two sequences using a whitelist. The actual data files will be reused from the input DB.
Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ rm -f test-1.filtered*
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" > filterlist.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> filterlist.txt
  $ ${BIN_DIR}/pancake dbfilter --block-size 0 --filter-list filterlist.txt --filter-type whitelist ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  V	0.1.0
  C	1
  F	0	test-1.seqdb.0.seq	1	1463	5852
  F	1	test-1.seqdb.1.seq	1	2996	11983
  F	2	test-1.seqdb.2.seq	1	6073	24292
  F	3	test-1.seqdb.3.seq	1	1277	5105
  F	4	test-1.seqdb.4.seq	1	4751	19001
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	1	0	2996	11983	1	0	11983
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	4	0	4751	19001	1	0	19001
  B	0	0	1	2996	11983
  B	1	1	2	4751	19001
  test-1.filtered.seqdb

Consolidate the DB.
Filter out two sequences using a whitelist. The new data files will be generated, and these should contain only the non-filtered sequences.
Also, unlike the input DB, here all the sequences will be stored in the same file.
  $ rm -f test-1.filtered*
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" > filterlist.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> filterlist.txt
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --filter-list filterlist.txt --filter-type whitelist --consolidate ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.4.seq > expected.seq
  > diff expected.seq test-1.filtered.seqdb.0.seq
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  V	0.1.0
  C	1
  F	0	test-1.filtered.seqdb.0.seq	2	7747	30984
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	0	2996	11983	1	0	11983
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	2996	4751	19001	1	0	19001
  B	0	0	1	2996	11983
  B	1	1	2	4751	19001
  test-1.filtered.seqdb
  test-1.filtered.seqdb.0.seq
##############################

####################
### Subsampling. ###
####################
Wrong subsampling type, should throw.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling something_wrong ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered 2>&1 | sed 's/^.*pancake //g'
  dbfilter ERROR: Unknown sampling type: 'something_wrong'.

No sampling, output should be the same as input.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling none ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  test-1.filtered.seqdb

Linear sampling and consolidation, normal valid case.
Pick the first ~45000bp from the input SeqDB, and create a new DB with consolidated sequences.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling linear --sample-bases 45000 --consolidate ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.0.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.3.seq > expected.seq
  > diff expected.seq test-1.filtered.seqdb.0.seq
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  V	0.1.0
  C	1
  F	0	test-1.filtered.seqdb.0.seq	4	11809	47232
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  B	0	0	1	1463	5852
  B	1	1	2	2996	11983
  B	2	2	3	6073	24292
  B	3	3	4	1277	5105
  test-1.filtered.seqdb
  test-1.filtered.seqdb.0.seq

Linear sampling and consolidation, with 0bp sampled bases. Output should have no sequences.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling linear --sample-bases 0 --consolidate ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  > wc -l test-1.filtered.seqdb.0.seq | awk '{ print $1 }'
  V	0.1.0
  C	1
  F	0	test-1.filtered.seqdb.0.seq	0	0	0
  test-1.filtered.seqdb
  test-1.filtered.seqdb.0.seq
  0

Random permutation and consolidation.
Filter out two sequences using a whitelist. The new data files will be generated, and these should contain only the non-filtered sequences.
Also, unlike the input DB, here all the sequences will be stored in the same file.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling random --sample-bases 25000 --random-seed 12345 --consolidate ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq > expected.seq
  > diff expected.seq test-1.filtered.seqdb.0.seq
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  V	0.1.0
  C	1
  F	0	test-1.filtered.seqdb.0.seq	2	9069	36275
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	0	2996	11983	1	0	11983
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	2996	6073	24292	1	0	24292
  B	0	0	1	2996	11983
  B	1	1	2	6073	24292
  test-1.filtered.seqdb
  test-1.filtered.seqdb.0.seq

Another random andom permutation and consolidation, but with a different seed.
Filter out two sequences using a whitelist. The new data files will be generated, and these should contain only the non-filtered sequences.
Also, unlike the input DB, here all the sequences will be stored in the same file.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling random --sample-bases 25000 --random-seed 7 --consolidate ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.0.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.3.seq > expected.seq
  > diff expected.seq test-1.filtered.seqdb.0.seq
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  V	0.1.0
  C	1
  F	0	test-1.filtered.seqdb.0.seq	3	8813	35249
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	1463	6073	24292	1	0	24292
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	7536	1277	5105	1	0	5105
  B	0	0	1	1463	5852
  B	1	1	2	6073	24292
  B	2	2	3	1277	5105
  test-1.filtered.seqdb
  test-1.filtered.seqdb.0.seq

Random sampling and consolidation, with 0bp sampled bases. Output should have no sequences.
  $ rm -f test-1.filtered*
  > ${BIN_DIR}/pancake dbfilter --block-size 0 --sampling random --sample-bases 0 --consolidate ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file test-1.filtered
  > cat test-1.filtered.seqdb
  > ls -1 test-1.filtered*
  > wc -l test-1.filtered.seqdb.0.seq | awk '{ print $1 }'
  V	0.1.0
  C	1
  F	0	test-1.filtered.seqdb.0.seq	0	0	0
  test-1.filtered.seqdb
  test-1.filtered.seqdb.0.seq
  0
