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
### Test tools for a SeqDB input without sequences, but with a header.
SeqDB test for empty input.
  $ rm -rf out; mkdir -p out
  > touch out/in.fasta
  > ${BIN_DIR}/pancake seqdb out/out out/in.fasta

SeedDB test for empty SeqDB, which still includes a header.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seeddb ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb out/out

Ovl-hifi test for empty SeqDB, which still includes a header.. This test is expected to throw, because the target block ID ("0") is explicitly specified, while there are no blocks in the input file.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake ovl-hifi ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: Invalid blockId (a). blockId = 0, blockLines.size() = 0

DBFilter test for empty SeqDB, which still includes a header..
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake dbfilter ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs out/out

Seqfetch test for empty SeqDB, which still includes a header..
  $ rm -rf out; mkdir -p out
  > touch out/in.txt
  > ${BIN_DIR}/pancake seqfetch out/out.fasta out/in.txt ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb

Seqdb-dump test for empty SeqDB, which still includes a header..
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb-dump ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb out/out.fasta

Seqdb-info test for empty SeqDB, which still includes a header..
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb-info ${PROJECT_DIR}/test-data/bugfixes/empty-seqdb/test-1-seqdb-with-no-seqs.seqdb --human | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	0.00	0	0.00	0.00	0.00	0.00	0.00	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0

Composite test for multiple tools, from an empty input FASTA.
  $ rm -rf out; mkdir -p out
  > touch out/in.fasta
  > ${BIN_DIR}/pancake seqdb out/out out/in.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  > ${BIN_DIR}/pancake dbfilter out/out out/out.filtered
  > touch out/to_fetch.txt
  > ${BIN_DIR}/pancake seqfetch out/out.fetched.fasta out/to_fetch.txt out/out.seqdb
  > wc -l out/out.fetched.fasta | awk '{ print $1 }'
  > ${BIN_DIR}/pancake seqdb-dump out/out.seqdb out/out.all.fasta
  > diff out/in.fasta out/out.all.fasta
  > ${BIN_DIR}/pancake seqdb-info out/out.seqdb --human | cut -f 2-
  ovl-hifi ERROR: Invalid blockId (a). blockId = 0, blockLines.size() = 0
  0
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	0.00	0	0.00	0.00	0.00	0.00	0.00	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0

### Test tools for a completely empty (zero-byte) SeqDB.
SeedDB test for a zero-byte SeqDB.
  $ rm -rf out; mkdir -p out
  > touch out/in.seqdb
  > ${BIN_DIR}/pancake seeddb out/in.seqdb out/out

Ovl-hifi test for a zero-byte SeqDB. This test is expected to throw, because the target block ID ("0") is explicitly specified, while there are no blocks in the input file.
  $ rm -rf out; mkdir -p out
  > touch out/in.seqdb
  > touch out/in.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/in out/in 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: Invalid blockId (a). blockId = 0, blockLines.size() = 0

DBFilter test for a zero-byte SeqDB.
  $ rm -rf out; mkdir -p out
  > touch out/in.seqdb
  > ${BIN_DIR}/pancake dbfilter out/in out/out

Seqfetch test for a zero-byte SeqDB.
  $ rm -rf out; mkdir -p out
  > touch out/in.seqdb
  > touch out/in.txt
  > ${BIN_DIR}/pancake seqfetch out/out.fasta out/in.txt out/in.seqdb

Seqdb-dump test for a zero-byte SeqDB.
  $ rm -rf out; mkdir -p out
  > touch out/in.seqdb
  > ${BIN_DIR}/pancake seqdb-dump out/in.seqdb out/out.fasta

Seqdb-info test for a zero-byte SeqDB.
  $ rm -rf out; mkdir -p out
  > touch out/in.seqdb
  > ${BIN_DIR}/pancake seqdb-info out/in.seqdb --human | cut -f 2-
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	0.00	0	0.00	0.00	0.00	0.00	0.00	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0

Composite test for multiple tools, from an empty input FASTA.
  $ rm -rf out; mkdir -p out
  > touch out/in.fasta
  > touch out/out.seqdb
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  > ${BIN_DIR}/pancake dbfilter out/out out/out.filtered
  > touch out/to_fetch.txt
  > ${BIN_DIR}/pancake seqfetch out/out.fetched.fasta out/to_fetch.txt out/out.seqdb
  > wc -l out/out.fetched.fasta | awk '{ print $1 }'
  > ${BIN_DIR}/pancake seqdb-dump out/out.seqdb out/out.all.fasta
  > diff out/in.fasta out/out.all.fasta
  > ${BIN_DIR}/pancake seqdb-info out/out.seqdb --human | cut -f 2-
  ovl-hifi ERROR: Invalid blockId (a). blockId = 0, blockLines.size() = 0
  0
  unit	total	num	min	max	avg	median	AUC	N10	N10_n	N25	N25_n	N50	N50_n	N75	N75_n	N90	N90_n	N100	N100_n
  bp	0.00	0	0.00	0.00	0.00	0.00	0.00	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0	0.00	0

### Test Pancake on combinations of SeqDBs where File, Sequence and/or Block lines are missing from the SeqDB file. In some of these cases,
### the tools should throw an ERROR, but not in all cases.
# F S B
# 1 1 1 FSB
# 1 1 0 FS
# 1 0 1 FB
# 1 0 0 F
# 0 1 1 SB
# 0 1 0 S
# 0 0 1 B
# 0 0 0 -
All three (File, Sequence and Block) lines are present in the input.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb > out/v1.FSB.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'

File and Sequence lines are present, but not Block lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[B]" > out/v2.FS.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'
  seeddb ERROR: There are no block specifications in the input SeqDB index file, but there are sequence lines listed.

File and Block lines are present, but not Sequence lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[S]" > out/v3.FB.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'
  seeddb ERROR: There are blocks specified in the input SeqDB index file, but there are no sequences listed.

File lines are present, but not Sequence or Block lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[SB]" > out/v4.F.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'

Sequence and Block lines are present, but not File lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[F]" > out/v5.SB.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'
  seeddb ERROR: There are no file specifications in the input SeqDB index file, but there are sequence lines listed.

Sequence lines are present, but not Block or File lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[FB]" > out/v6.S.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'
  seeddb ERROR: There are no file specifications in the input SeqDB index file, but there are sequence lines listed.

Block lines are present, but not Sequence or File lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[FS]" > out/v7.B.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'
  seeddb ERROR: There are no file specifications in the input SeqDB index file, but there are block lines listed.

There are no File, Sequence or Block lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > cat out/out.seqdb | grep -v "^[FSB]" > out/v8.no_FSB.seqdb
  > ${BIN_DIR}/pancake seeddb out/v*.seqdb out/out 2>&1 | sed 's/.*pancake //g'

### Test Pancake on combinations of SeedDBs where File, Seeds and/or Block lines are missing from the SeedDB file. In some of these cases,
### the tool should throw an ERROR, but not in all cases.
# F S B
# 1 1 1 FSB
# 1 1 0 FS
# 1 0 1 FB
# 1 0 0 F
# 0 1 1 SB
# 0 1 0 S
# 0 0 1 B
# 0 0 0 -
All three (File, Seed and Block) lines are present in the input.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'

File and Seed lines are present, but not Block lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[B]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: There are no block specifications in the input SeedDB index file, but there are seed lines listed.

File and Block lines are present, but not Seed lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[S]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: There are blocks specified in the input SeedDB index file, but there are no seed lines listed.

File lines are present, but not Seed or Block lines.
In this case, the SeedDB itself is not formatted wrongly, but it does not match the accompanying SeqDB (i.e.
it is missing the target block required for overlapping).
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[SB]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: Invalid blockId (SeedDBReader). blockId = 0, blocks.size() = 0

Seed and Block lines are present, but not File lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[F]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: There are no file specifications in the input SeedDB index file, but there are seed lines listed.

Seed lines are present, but not Block or File lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[FB]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: There are no file specifications in the input SeedDB index file, but there are seed lines listed.

Block lines are present, but not Seed or File lines.
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[FS]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: There are no file specifications in the input SeedDB index file, but there are block lines listed.

There are no File, Seed or Block lines.
In this case, the SeedDB itself is not formatted wrongly, but it does not match the accompanying SeqDB (i.e.
it is missing the target block required for overlapping).
  $ rm -rf out; mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/out ${PROJECT_DIR}/test-data/seqdb-writer/in-2-small.fasta
  > ${BIN_DIR}/pancake seeddb out/out.seqdb out/out.vanilla 2>&1 | sed 's/.*pancake //g'
  > cat out/out.vanilla.seeddb | grep -v "^[FSB]" > out/out.seeddb
  > ${BIN_DIR}/pancake ovl-hifi out/out out/out 0 0 0 2>&1 | sed 's/.*pancake //g'
  ovl-hifi ERROR: Invalid blockId (SeedDBReader). blockId = 0, blocks.size() = 0
############################
