Test fetching of sequences from two input files, one FASTA and one FASTQ.
Output in FASTA.
  $ rm -f out.* test.in.*
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta > test.in.1.fasta
  > tail -n 8 ${TEST_DATA_DIR}/seqdb-writer/in.fastq > test.in.2.fastq
  > samtools faidx test.in.1.fasta
  > samtools fqidx test.in.2.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt test.in.1.fasta test.in.2.fastq
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta | tail -n 2 > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from the SeqDB.
Output in FASTA.
  $ rm -f out.* test.in.*
  > ${BIN_DIR}/pancake seqdb test.in.1 ${TEST_DATA_DIR}/seqdb-writer/in.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt test.in.1.seqdb
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta | tail -n 2 > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from a BAM file.
Output in FASTA.
  $ rm -f out.* test.in.*
  > echo "ref1/3/0_42176" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt ${TEST_DATA_DIR}/seqdb-writer/bam/subreads1.bam
  > samtools fasta ${TEST_DATA_DIR}/seqdb-writer/bam/subreads1.bam | head -n 6 | tail -n 2 > expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from two input files, one FASTA and one FASTQ.
Output in FASTA.
  $ rm -f out.* test.in.*
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta > test.in.1.fasta
  > tail -n 8 ${TEST_DATA_DIR}/seqdb-writer/in.fastq > test.in.2.fastq
  > samtools faidx test.in.1.fasta
  > samtools fqidx test.in.2.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852" > test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt test.in.1.fasta test.in.2.fastq
  > grep ">" ${TEST_DATA_DIR}/seqdb-writer/in.fasta | sort > expected.headers
  > grep ">" out.fasta | sort > results.headers
  > diff expected.headers results.headers | wc -l | awk '{ print $1 }'
  0

Output in FASTQ.
  $ rm -f out.* test.in.*
  > cp ${TEST_DATA_DIR}/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fastq out.fastq test.in.list.txt test.in.1.fastq
  > head -n 12 ${TEST_DATA_DIR}/seqdb-writer/in.fastq | tail -n 4 > expected.fastq
  > diff expected.fastq out.fastq | wc -l | awk '{ print $1 }'
  0

Test fetching of sequences from a FOFN.
  $ rm -f out.* test.in.*
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta > test.in.1.fasta
  > tail -n 8 ${TEST_DATA_DIR}/seqdb-writer/in.fastq > test.in.2.fastq
  > samtools faidx test.in.1.fasta
  > samtools fqidx test.in.2.fastq
  > echo "test.in.1.fasta" > input.fofn
  > echo "test.in.2.fastq" >> input.fofn
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3005/0_5852" > test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3414/0_11983" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105" >> test.in.list.txt
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta out.fasta test.in.list.txt input.fofn
  > grep ">" ${TEST_DATA_DIR}/seqdb-writer/in.fasta | sort > expected.headers
  > grep ">" out.fasta | sort > results.headers
  > diff expected.headers results.headers | wc -l | awk '{ print $1 }'
  0

Test fetching using alias. Output in FASTA.
  $ rm -f out.* test.in.*
  > cp ${TEST_DATA_DIR}/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > ${BIN_DIR}/pancake seqdb test.in.1 test.in.1.fastq
  > echo "3" > test.in.list.txt
  > echo "4" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta --alias test.in.1.seqdb out.fasta test.in.list.txt test.in.1.fastq
  > echo ">m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3981/0_5105" > expected.headers
  > echo ">m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/4028/0_19001" >> expected.headers
  > cat expected.headers | sort > expected.headers.sorted
  > grep ">" out.fasta | sort > results.headers
  > diff expected.headers.sorted results.headers | wc -l | awk '{ print $1 }'
  0

Nonexistent sequence, but do not fail.
  $ rm -f out.* test.in.*
  > cp ${TEST_DATA_DIR}/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "some_read_which_does_not_exist" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fastq out.fastq test.in.list.txt test.in.1.fastq
  > wc -l out.fastq | awk '{ print $1 }'
  0

Nonexistent sequence, FAIL when it cannot be found.
  $ rm -f out.* test.in.*
  > cp ${TEST_DATA_DIR}/seqdb-writer/in.fastq test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "some_read_which_does_not_exist" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --fail --out-fmt fastq out.fastq test.in.list.txt test.in.1.fastq 2>&1 | sed 's/.*pancake //g'
  seqfetch ERROR: Not all queries were found in the provided input files! Found sequences: 0 / 1 .

Test writing the sequence ID instead of the header.
Output in FASTA.
  $ rm -f out.* test.in.*
  > ${BIN_DIR}/pancake seqdb test.in.1 ${TEST_DATA_DIR}/seqdb-writer/in.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --out-fmt fasta --write-ids out.fasta test.in.list.txt test.in.1.seqdb
  > echo ">000000002" > expected.fasta
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta | tail -n 1 >> expected.fasta
  > diff expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  0

Test that the --write-ids feature should fail if the input is not SeqDB.
The list of input files consists of one SeqDB file and one FASTA file. The FASTA wille should throw.
  $ rm -f out.* test.in.*
  > head -n 6 ${TEST_DATA_DIR}/seqdb-writer/in.fasta > test.in.1.fasta
  > samtools faidx test.in.1.fasta
  > ${BIN_DIR}/pancake seqdb test.in.2 ${TEST_DATA_DIR}/seqdb-writer/in.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --write-ids --out-fmt fasta out.fasta test.in.list.txt test.in.2.seqdb test.in.1.fasta 2>&1 | sed 's/.*pancake //g'
  seqfetch ERROR: Cannot use the --write-ids option with input files which are not in the SeqDB format. Offending file: 'test.in.1.fasta'.

Test that the --write-ids feature should fail if the input is not SeqDB.
Here, the offending input is in FASTQ.
  $ rm -f out.* test.in.*
  > tail -n 8 ${TEST_DATA_DIR}/seqdb-writer/in.fastq > test.in.1.fastq
  > samtools fqidx test.in.1.fastq
  > echo "m141013_011508_sherri_c100709962550000001823135904221533_s1_p0/3820/0_24292" > test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --write-ids --out-fmt fasta out.fasta test.in.list.txt test.in.1.fastq 2>&1 | sed 's/.*pancake //g'
  seqfetch ERROR: Cannot use the --write-ids option with input files which are not in the SeqDB format. Offending file: 'test.in.1.fastq'.

####################
### Testing HPC. ###
####################

Test fetching with homopolymer compression.
Output in FASTA.
  $ rm -f out.* test.in.*
  > echo "seq-1-single_base" > test.in.list.txt
  > echo "seq-2-single_polyA_5x" >> test.in.list.txt
  > echo "seq-3-simple_seq" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --use-hpc --out-fmt fasta out.fasta test.in.list.txt ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-expected.fasta out.fasta
  0

Test the '--use-rle' option. The fetched sequences should be identical to the original ones because we are not
specifying the '--use-hpc' here, but there should be an additional file written that contains the run length encoding, out.fasta.rle.
Output in FASTA.
  $ rm -f out.* test.in.*
  > echo "seq-1-single_base" > test.in.list.txt
  > echo "seq-2-single_polyA_5x" >> test.in.list.txt
  > echo "seq-3-simple_seq" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --use-rle --out-fmt fasta out.fasta test.in.list.txt ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta out.fasta | wc -l | awk '{ print $1 }'
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-expected.rle out.fasta.rle | wc -l | awk '{ print $1 }'
  0
  0

Test the '--use-rle' and '--use-hpc' options together. The fetched sequences should be homopolymer compressed, and
there should be an additional file written that contains the run length encoding, out.fasta.rle.
Output in FASTA.
  $ rm -f out.* test.in.*
  > echo "seq-1-single_base" > test.in.list.txt
  > echo "seq-2-single_polyA_5x" >> test.in.list.txt
  > echo "seq-3-simple_seq" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --use-hpc --use-rle --out-fmt fasta out.fasta test.in.list.txt ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-expected.fasta out.fasta | wc -l | awk '{ print $1 }'
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-expected.rle out.fasta.rle | wc -l | awk '{ print $1 }'
  0
  0

Test fetching with homopolymer compression and run-length encoding, but try outputting to stdout.
This should fail.
  $ rm -f out.* test.in.*
  > echo "seq-1-single_base" > test.in.list.txt
  > echo "seq-2-single_polyA_5x" >> test.in.list.txt
  > echo "seq-3-simple_seq" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --use-rle --out-fmt fasta - test.in.list.txt ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta 2>&1 | sed 's/.*pancake //g'
  seqfetch ERROR: Cannot output to sequences to stdout and write a .rle file. Please specify a concrete output file.

Test the '--use-rle' option without '--use-hpc' and output to FASTQ. This should be valid and the output sequences
the same as input.
Output in FASTA.
  $ rm -f out.* test.in.*
  > echo "seq-1-single_base" > test.in.list.txt
  > echo "seq-2-single_polyA_5x" >> test.in.list.txt
  > echo "seq-3-simple_seq" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --use-rle --out-fmt fastq out.fastq test.in.list.txt ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta
  > diff ${TEST_DATA_DIR}/seqfetch/test-1-hpc-expected.rle out.fastq.rle | wc -l | awk '{ print $1 }'
  > wc -l out.fastq | awk '{ print $1 }'
  0
  12

Test fetching with homopolymer compression and run-length encoding, but try outputting to a FASTQ file.
This should fail.
  $ rm -f out.* test.in.*
  > echo "seq-1-single_base" > test.in.list.txt
  > echo "seq-2-single_polyA_5x" >> test.in.list.txt
  > echo "seq-3-simple_seq" >> test.in.list.txt
  > ${BIN_DIR}/pancake seqfetch --use-hpc --out-fmt fastq out.fastq test.in.list.txt ${TEST_DATA_DIR}/seqfetch/test-1-hpc-in.fasta 2>&1 | sed 's/.*pancake //g'
  seqfetch ERROR: Fastq output format is not supported with the homopolymer compression option.