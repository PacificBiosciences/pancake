Test the seqdb-info tool on a small sample DB.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/seqdb-info/test-1-compressed-each-seq-one-block-and-file.seqdb --human | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	66233.00	5	5105.00	24292.00	13246.60	11983.00	17439.02	24292.00	1	24292.00	1	19001.00	2	11983.00	3	5852.00	4	5105.00	5

Change the output unit to kilobases.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/seqdb-info/test-1-compressed-each-seq-one-block-and-file.seqdb --human --unit kbp | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  kbp	66.23	5	5.11	24.29	13.25	11.98	17.44	24.29	1	24.29	1	19.00	2	11.98	3	5.85	4	5.11	5

Change the output unit to megabases.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/seqdb-info/test-1-compressed-each-seq-one-block-and-file.seqdb --human --unit Mbp | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  Mbp	0.07	5	0.01	0.02	0.01	0.01	0.02	0.02	1	0.02	1	0.02	2	0.01	3	0.01	4	0.01	5

Change the output unit to gigabases.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/seqdb-info/test-1-compressed-each-seq-one-block-and-file.seqdb --human --unit Gbp | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  Gbp	0.00	5	0.00	0.00	0.00	0.00	0.00	0.00	1	0.00	1	0.00	2	0.00	3	0.00	4	0.00	5

Write the stats as a JSON output.
  $ rm -f out.*
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/seqdb-info/test-1-compressed-each-seq-one-block-and-file.seqdb > out.json
  > diff ${PROJECT_DIR}/test-data/seqdb-info/test-1-compressed-each-seq-one-block-and-file.out.json out.json
