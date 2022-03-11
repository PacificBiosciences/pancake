Test construction of the DB from a small test FASTA file. Store each sequence into a separate 2-bit compressed file.
  $ rm -f test-1.seqdb*
  > ${BIN_DIR}/pancake seqdb test-1 ${TEST_DATA_DIR}/seqdb-writer/in.fasta --block-size 0 --buffer-size 1024 --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.0.seq test-1.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq test-1.seqdb.1.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq test-1.seqdb.2.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.3.seq test-1.seqdb.3.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.4.seq test-1.seqdb.4.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.seqdb
  > ls -1 test-1.seqdb*
  test-1.seqdb
  test-1.seqdb.0.seq
  test-1.seqdb.1.seq
  test-1.seqdb.2.seq
  test-1.seqdb.3.seq
  test-1.seqdb.4.seq

Test construction of the DB from a small test FASTA file, and split blocks but this time more than 1 sequence per block.
  $ rm -f test-1.seqdb*
  > ${BIN_DIR}/pancake seqdb test-11-split-blocks ${TEST_DATA_DIR}/seqdb-writer/in.fasta --block-size 0.015 --buffer-size 1024 --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-11-split-blocks.seqdb test-11-split-blocks.seqdb
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-11-split-blocks.seqdb.0.seq test-11-split-blocks.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-11-split-blocks.seqdb.1.seq test-11-split-blocks.seqdb.1.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-11-split-blocks.seqdb.2.seq test-11-split-blocks.seqdb.2.seq
  > ls -1 test-11-split-blocks.seqdb*
  test-11-split-blocks.seqdb
  test-11-split-blocks.seqdb.0.seq
  test-11-split-blocks.seqdb.1.seq
  test-11-split-blocks.seqdb.2.seq

FASTQ input.
Test construction of the DB from a small test FASTQ file. Store each sequence into a separate 2-bit compressed file.
This test is exactly the same as the previous one, the only difference is that the input is in the FASTQ format instead of FASTA.
That is why we can reuse the output files from the previous test, even though the input is now "in.fastq" instead of "in.fasta".
Note: the FASTQ here contains identical sequences to the "in.fasta", and the quality values are just dummy values added to create a FASTQ.
  $ rm -f test-1.seqdb*
  > ${BIN_DIR}/pancake seqdb test-1 ${TEST_DATA_DIR}/seqdb-writer/in.fastq --block-size 0 --buffer-size 1024 --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.0.seq test-1.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq test-1.seqdb.1.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq test-1.seqdb.2.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.3.seq test-1.seqdb.3.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.4.seq test-1.seqdb.4.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.seqdb
  > ls -1 test-1.seqdb*
  test-1.seqdb
  test-1.seqdb.0.seq
  test-1.seqdb.1.seq
  test-1.seqdb.2.seq
  test-1.seqdb.3.seq
  test-1.seqdb.4.seq

FOFN input.
Test construction of the DB from a small test FASTQ file. Store each sequence into a separate 2-bit compressed file.
This test is exactly the same as the previous one, the only difference is that the input is in the FOFN format instead of FASTA.
That is why we can reuse the output files from the previous test, even though the input is now "in.fofn" instead of "in.fasta".
  $ rm -f test-1.seqdb*
  > echo "${TEST_DATA_DIR}/seqdb-writer/in.fasta" > in.fofn
  > ${BIN_DIR}/pancake seqdb test-1 in.fofn --block-size 0 --buffer-size 1024 --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.0.seq test-1.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq test-1.seqdb.1.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq test-1.seqdb.2.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.3.seq test-1.seqdb.3.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.4.seq test-1.seqdb.4.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.seqdb
  > ls -1 test-1.seqdb*
  test-1.seqdb
  test-1.seqdb.0.seq
  test-1.seqdb.1.seq
  test-1.seqdb.2.seq
  test-1.seqdb.3.seq
  test-1.seqdb.4.seq

Same as before, but test writing to a different folder. The file paths should be local.
  $ mkdir -p out
  > ${BIN_DIR}/pancake seqdb out/test-1 ${TEST_DATA_DIR}/seqdb-writer/in.fasta --block-size 0 --buffer-size 1024 --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.0.seq out/test-1.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.1.seq out/test-1.seqdb.1.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.2.seq out/test-1.seqdb.2.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.3.seq out/test-1.seqdb.3.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.4.seq out/test-1.seqdb.4.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-1-compressed-each-seq-one-block-and-file.seqdb test-1.seqdb
  > ls -1 out/test-1.seqdb*
  out/test-1.seqdb
  out/test-1.seqdb.0.seq
  out/test-1.seqdb.1.seq
  out/test-1.seqdb.2.seq
  out/test-1.seqdb.3.seq
  out/test-1.seqdb.4.seq

Write the SeqDB like before, but store all sequences into the same single 2-bit compressed file.
  $ rm -f out.seqdb*
  > ${BIN_DIR}/pancake seqdb out ${TEST_DATA_DIR}/seqdb-writer/in.fasta --block-size 1024 --buffer-size 1024
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.*.seq > expected.seq
  > cat out.seqdb
  > diff expected.seq out.seqdb.0.seq
  > ls -1 out.seqdb*
  V	0.1.0
  C	1
  F	0	out.seqdb.0.seq	5	16560	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
  B	0	0	5	16560	66233
  out.seqdb
  out.seqdb.0.seq

Using any size buffer should not impact the results, only perhaps the speed of writing. This will write out one read at a time.
  $ rm -f out.seqdb*
  > ${BIN_DIR}/pancake seqdb out ${TEST_DATA_DIR}/seqdb-writer/in.fasta --block-size 1024 --buffer-size 0
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-1.seqdb.*.seq > expected.seq
  > cat out.seqdb
  > diff expected.seq out.seqdb.0.seq
  > ls -1 out.seqdb*
  V	0.1.0
  C	1
  F	0	out.seqdb.0.seq	5	16560	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	1463	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	1463	2996	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	4459	6073	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	10532	1277	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	11809	4751	19001	1	0	19001
  B	0	0	5	16560	66233
  out.seqdb
  out.seqdb.0.seq

Uncompressed SeqDB construction, each sequence in one block.
  $ rm -f test-3.seqdb*
  > ${BIN_DIR}/pancake seqdb test-3 ${TEST_DATA_DIR}/seqdb-writer/in.fasta --compression 0 --block-size 0 --buffer-size 1024 --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-3.seqdb.0.seq test-3.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-3.seqdb.1.seq test-3.seqdb.1.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-3.seqdb.2.seq test-3.seqdb.2.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-3.seqdb.3.seq test-3.seqdb.3.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-3.seqdb.4.seq test-3.seqdb.4.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-3-uncompressed-each-seq-one-block-and-file.seqdb test-3.seqdb
  > ls -1 test-3.seqdb*
  test-3.seqdb
  test-3.seqdb.0.seq
  test-3.seqdb.1.seq
  test-3.seqdb.2.seq
  test-3.seqdb.3.seq
  test-3.seqdb.4.seq

Uncompressed SeqDB construction, single block.
  $ rm -f out.seqdb*
  > ${BIN_DIR}/pancake seqdb out ${TEST_DATA_DIR}/seqdb-writer/in.fasta --compression 0 --block-size 1024 --buffer-size 1024
  > cat ${TEST_DATA_DIR}/seqdb-writer/test-3.seqdb.*.seq > expected.seq
  > cat out.seqdb
  > diff expected.seq out.seqdb.0.seq
  > ls -1 out.seqdb*
  V	0.1.0
  C	0
  F	0	out.seqdb.0.seq	5	66233	66233
  S	0	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852	0	0	5852	5852	1	0	5852
  S	1	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983	0	5852	11983	11983	1	0	11983
  S	2	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292	0	17835	24292	24292	1	0	24292
  S	3	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105	0	42127	5105	5105	1	0	5105
  S	4	m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001	0	47232	19001	19001	1	0	19001
  B	0	0	5	66233	66233
  out.seqdb
  out.seqdb.0.seq

Create only one sequence file because all sequences fit into one block.
  $ rm -f test-6.seqdb*
  > ${BIN_DIR}/pancake seqdb test-6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta --split-blocks
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-6.seqdb.0.seq test-6.seqdb.0.seq
  > diff ${TEST_DATA_DIR}/seqdb-writer/test-6.seqdb test-6.seqdb
  > ls -1 test-6.seqdb*
  test-6.seqdb
  test-6.seqdb.0.seq

Test an unsupported input format.
  $ rm -f test-1.seqdb*
  > ${BIN_DIR}/pancake seqdb test-1 ${TEST_DATA_DIR}/seqdb-writer/README.md --block-size 0 --buffer-size 1024 --split-blocks 2>&1 | grep "Unknown input file extension for file" | wc -l | awk '{ print $1 }'
  1

Test BAM input.
  $ rm -f test-input-bam.seqdb*
  > ${BIN_DIR}/pancake seqdb test-input-bam ${TEST_DATA_DIR}/seqdb-writer/bam/subreads1.bam --block-size 1024 --buffer-size 1024 --split-blocks
  > cat test-input-bam.seqdb
  > ls -1 test-input-bam.seqdb*
  V	0.1.0
  C	1
  F	0	test-input-bam.seqdb.0.seq	5	52727	210907
  S	0	ref1/1/0_42148	0	0	10537	42148	1	0	42148
  S	1	ref1/2/0_42124	0	10537	10531	42124	1	0	42124
  S	2	ref1/3/0_42176	0	21068	10544	42176	1	0	42176
  S	3	ref1/4/0_42179	0	31612	10545	42179	1	0	42179
  S	4	ref1/5/0_42280	0	42157	10570	42280	1	0	42280
  B	0	0	5	52727	210907
  test-input-bam.seqdb
  test-input-bam.seqdb.0.seq

Test XML input.
  $ rm -f test-input-xml.seqdb*
  > ${BIN_DIR}/pancake seqdb test-input-xml ${TEST_DATA_DIR}/seqdb-writer/bam/subreadset.xml --block-size 1024 --buffer-size 1024 --split-blocks
  > cat test-input-xml.seqdb
  > ls -1 test-input-xml.seqdb*
  V	0.1.0
  C	1
  F	0	test-input-xml.seqdb.0.seq	15	152664	610647
  S	0	ref1/1/0_42148	0	0	10537	42148	1	0	42148
  S	1	ref1/2/0_42124	0	10537	10531	42124	1	0	42124
  S	2	ref1/3/0_42176	0	21068	10544	42176	1	0	42176
  S	3	ref1/4/0_42179	0	31612	10545	42179	1	0	42179
  S	4	ref1/5/0_42280	0	42157	10570	42280	1	0	42280
  S	5	ref2/1/0_42148	0	52727	10537	42148	1	0	42148
  S	6	ref2/2/0_42124	0	63264	10531	42124	1	0	42124
  S	7	ref2/3/0_42176	0	73795	10544	42176	1	0	42176
  S	8	ref2/4/0_42179	0	84339	10545	42179	1	0	42179
  S	9	ref2/5/0_42280	0	94884	10570	42280	1	0	42280
  S	10	ref1/158/0_23696	0	105454	5924	23696	1	0	23696
  S	11	ref1/159/0_42173	0	111378	10544	42173	1	0	42173
  S	12	ref1/160/0_42090	0	121922	10523	42090	1	0	42090
  S	13	ref1/161/0_38694	0	132445	9674	38694	1	0	38694
  S	14	ref1/162/0_42180	0	142119	10545	42180	1	0	42180
  B	0	0	15	152664	610647
  test-input-xml.seqdb
  test-input-xml.seqdb.0.seq

Test filtered XML input.
  $ rm -f test-input-xml.seqdb*
  > ${BIN_DIR}/pancake seqdb test-input-xml ${TEST_DATA_DIR}/seqdb-writer/bam/filtered.subreadset.xml --block-size 1024 --buffer-size 1024 --split-blocks
  > cat test-input-xml.seqdb
  > ls -1 test-input-xml.seqdb*
  V	0.1.0
  C	1
  F	0	test-input-xml.seqdb.0.seq	11	116004	464009
  S	0	ref1/1/0_42148	0	0	10537	42148	1	0	42148
  S	1	ref1/3/0_42176	0	10537	10544	42176	1	0	42176
  S	2	ref1/4/0_42179	0	21081	10545	42179	1	0	42179
  S	3	ref1/5/0_42280	0	31626	10570	42280	1	0	42280
  S	4	ref2/1/0_42148	0	42196	10537	42148	1	0	42148
  S	5	ref2/3/0_42176	0	52733	10544	42176	1	0	42176
  S	6	ref2/4/0_42179	0	63277	10545	42179	1	0	42179
  S	7	ref2/5/0_42280	0	73822	10570	42280	1	0	42280
  S	8	ref1/159/0_42173	0	84392	10544	42173	1	0	42173
  S	9	ref1/160/0_42090	0	94936	10523	42090	1	0	42090
  S	10	ref1/162/0_42180	0	105459	10545	42180	1	0	42180
  B	0	0	11	116004	464009
  test-input-xml.seqdb
  test-input-xml.seqdb.0.seq
