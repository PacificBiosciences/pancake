########################################################
### Test mapping of reads to reference with flagging ###
### of secondary alignments.                         ###
########################################################
The read maps to ctg.000000F fully, and only partially to ctg.000000F. The partial mapping is too short and
is not flagged as secondary.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/mapping/test-1-no-secondary-aln.reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake seqdb ref ${TEST_DATA_DIR}/mapping/test-1-no-secondary-aln.ref.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 ref.seqdb ref
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 ref reads 0 0 0 --out-fmt paf --mark-secondary
  hap1/186/0_12695	12695	0	12695	+	ctg.000000F:114000-128000	14001	602	13301	12695	12699	60	tp:A:P	NM:i:25	IT:f:99.8031	SC:i:-12670	AT:Z:c	BT:Z:C	VQ:Z:*	VT:Z:*

The read maps to ctg.000000F fully, and only partially to ctg.000000F. The partial mapping is too short and
is not flagged as secondary.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/mapping/test-2-secondary-aln.reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake seqdb ref ${TEST_DATA_DIR}/mapping/test-2-secondary-aln.ref.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 ref.seqdb ref
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 ref reads 0 0 0 --out-fmt paf --mark-secondary
  hap1/1/0_11225	11225	0	11225	-	ctg.000000F:124200-136500	12301	697	11922	11225	11225	60	tp:A:P	NM:i:18	IT:f:99.8396	SC:i:-11207	AT:Z:c	BT:Z:C	VQ:Z:*	VT:Z:*
  hap1/1/0_11225	11225	0	11225	+	ctg.000001F:66000-77700	11701	213	11435	11225	11222	60	tp:A:S	NM:i:64	IT:f:99.4298	SC:i:-11158	AT:Z:c	BT:Z:C	VQ:Z:*	VT:Z:*

#######################################
### Test trimming of CIGAR strings. ###
#######################################
Simple case where no trimming should be applied. The trailing "1D" is allowed because the option "--trim-to-first-match" is not specified.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --no-indels --use-hpc --min-map-len 0 --min-anchor-span 100 --write-cigar --trim --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 100.00 0 0 5145 6298 0 496 5640 5640 5 221=1D649=1D296=1I169=1I167=1I3640=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 100.00 0 0 5026 6298 0 1470 6497 6497 5 5026=1D
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3 5008=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5 1288=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3 952=

Trim to first match
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --no-indels --use-hpc --min-map-len 0 --min-anchor-span 100 --write-cigar --trim --trim-to-first-match --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 100.00 0 0 5145 6298 0 496 5640 5640 5 221=1D649=1D296=1I169=1I167=1I3640=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 100.00 0 0 5026 6298 0 1470 6496 6497 u 5026=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3 5008=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5 1288=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3 952=
