#################################################
### Duplicate alignments reported.            ###
#################################################
This pair tended to generate multiple copies of the same overlap. It used to look like this:
# mo_m64001_190914_015449/117639449/ccs fa_m64001_190914_015449/100008058/ccs -16220 98.4717 0 2184 18738 18738 0 0 16473 18083 3 5 u * * * *
# mo_m64001_190914_015449/117639449/ccs fa_m64001_190914_015449/100008058/ccs -16220 98.4717 0 2184 18738 18738 0 0 16473 18083 3 5 u * * * *
# fa_m64001_190914_015449/100008058/ccs mo_m64001_190914_015449/117639449/ccs -16220 98.4717 0 0 16473 18083 0 2184 18738 18738 5 3 u * * * *
# fa_m64001_190914_015449/100008058/ccs mo_m64001_190914_015449/117639449/ccs -16220 98.4717 0 0 16473 18083 0 2184 18738 18738 5 3 u * * * *
  $ ${BIN_DIR}/pancake seqdb --compression 1 reads ${PROJECT_DIR}/test-data/hifi-ovl/robustness/reads.robustness.duplicate.2.single_pair.fasta
  > ${BIN_DIR}/pancake seeddb -k 28 -w 80 --space 1 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --skip-sym --write-rev --out-fmt ipa
  mo_m64001_190914_015449/117639449/ccs fa_m64001_190914_015449/100008058/ccs -16220 98.4717 0 2184 18738 18738 0 0 16473 18083 3 5 u * * * *
  fa_m64001_190914_015449/100008058/ccs mo_m64001_190914_015449/117639449/ccs -16220 98.4717 0 0 16473 18083 0 2184 18738 18738 5 3 u * * * *

Another test case for duplicates.
  $ ${BIN_DIR}/pancake seqdb --compression 1 reads ${PROJECT_DIR}/test-data/hifi-ovl/robustness/reads.robustness.duplicate.5.single_pair.dup_with_different_idt.fasta
  > ${BIN_DIR}/pancake seeddb -k 28 -w 80 --space 1 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --min-idt 96.0 --skip-sym --write-rev --out-fmt ipa
  mo_m64001_190914_015449/49284690/ccs fa_m64001_190914_015449/100010148/ccs -15386 97.9491 0 0 15798 24551 1 0 15710 18018 5 5 u * * * *
  fa_m64001_190914_015449/100010148/ccs mo_m64001_190914_015449/49284690/ccs -15386 97.9491 0 0 15710 18018 1 0 15798 24551 5 5 u * * * *

This dataset used to generate duplicated overlaps, but with slightly different identities, even though the coordinates were the same.
The output used to look like this:
# mo_m64001_190914_015449/14550518/ccs fa_m64001_190914_015449/100336406/ccs -10558 96.7229 0 0 11016 19863 0 6562 17481 17481 5 3 u * * * *
# fa_m64001_190914_015449/100336406/ccs mo_m64001_190914_015449/14550518/ccs -10558 96.7229 0 6562 17481 17481 0 0 11016 19863 3 5 u * * * *
# mo_m64001_190914_015449/14550518/ccs fa_m64001_190914_015449/100336406/ccs -10528 96.4506 0 0 11016 19863 0 6562 17481 17481 5 3 u * * * *
# fa_m64001_190914_015449/100336406/ccs mo_m64001_190914_015449/14550518/ccs -10528 96.4506 0 6562 17481 17481 0 0 11016 19863 3 5 u * * * *
  $ ${BIN_DIR}/pancake seqdb --compression 1 reads ${PROJECT_DIR}/test-data/hifi-ovl/robustness/reads.robustness.duplicate.9.single_pair.dup_with_different_idt.fasta
  > ${BIN_DIR}/pancake seeddb -k 28 -w 80 --space 1 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --min-idt 96.0 --skip-sym --write-rev --out-fmt ipa
  mo_m64001_190914_015449/14550518/ccs fa_m64001_190914_015449/100336406/ccs -10558 96.7229 0 0 11016 19863 0 6562 17481 17481 5 3 u * * * *
  fa_m64001_190914_015449/100336406/ccs mo_m64001_190914_015449/14550518/ccs -10558 96.7229 0 6562 17481 17481 0 0 11016 19863 3 5 u * * * *

#################################################
### Alignment loss due to spurious seed hits. ###
#################################################
This pair of reads presented a problem for mapping/overlapping, where one of the pair combinations had a bad seed hit,
which in turn cause the alignment to be wonky and filtered out in one strand.
Before the LIS was implemented, the following test would produce 0 alignments.
  $ ${BIN_DIR}/pancake seqdb --compression 1 ref ${PROJECT_DIR}/test-data/hifi-ovl/robustness/reads.robustness.spurious_hits.1.ref.fasta
  > ${BIN_DIR}/pancake seeddb -k 28 -w 80 --space 1 ref.seqdb ref
  > ${BIN_DIR}/pancake seqdb --compression 1 reads ${PROJECT_DIR}/test-data/hifi-ovl/robustness/reads.robustness.spurious_hits.1.reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 28 -w 80 --space 1 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --out-fmt ipa --one-hit-per-target ref reads 0 0 0 --out-fmt ipa --mark-secondary
  m64030_190330_071939/109841220/ccs m64030_190330_071939/94044350/ccs -9386 99.8619 0 0 9399 11041 1 0 9411 11621 5 5 u * * * *
