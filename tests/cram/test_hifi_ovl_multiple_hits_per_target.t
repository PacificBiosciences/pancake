#############################################
### Compare different ways of filtering   ###
### multiple hits per query/target pair.  ###
#############################################
All hits per query/target pair.
  $ rm -rf out && mkdir -p out
  > in_fasta=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.in.reads.plasmid_5kbp.ecoli_w_bc1087.fasta
  > expected_ovl=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.out.all_hits.ovl
  > ${BIN_DIR}/pancake seqdb out/reads ${in_fasta}
  > ${BIN_DIR}/pancake seeddb out/reads.seqdb out/reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --num-threads 1 --skip-sym --write-rev --out-fmt ipa --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 out/reads out/reads 0 0 0 --out-fmt m4 | sort -k 1,1 > out/out.ovl
  > diff ${expected_ovl} out/out.ovl

Marking secondary alignments and keeping only primary. Not well suited for overlapping because all overlaps are secondary, except
potential tandem repeats and circular matches.
  $ rm -rf out && mkdir -p out
  > in_fasta=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.in.reads.plasmid_5kbp.ecoli_w_bc1087.fasta
  > expected_ovl=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.out.mark_secondary.ovl
  > ${BIN_DIR}/pancake seqdb out/reads ${in_fasta}
  > ${BIN_DIR}/pancake seeddb out/reads.seqdb out/reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --num-threads 1 --skip-sym --write-rev --out-fmt ipa --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 out/reads out/reads 0 0 0 --out-fmt m4 --mark-secondary --bestn 0 | sort -k 1,1 > out/out.ovl
  > diff ${expected_ovl} out/out.ovl

Classic one-hit-per-target. This will not capture the circular matches.
  $ rm -rf out && mkdir -p out
  > in_fasta=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.in.reads.plasmid_5kbp.ecoli_w_bc1087.fasta
  > expected_ovl=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.out.one_hit_per_target.ovl
  > ${BIN_DIR}/pancake seqdb out/reads ${in_fasta}
  > ${BIN_DIR}/pancake seeddb out/reads.seqdb out/reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --num-threads 1 --skip-sym --write-rev --out-fmt ipa --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 out/reads out/reads 0 0 0 --out-fmt m4 --one-hit-per-target | sort -k 1,1 > out/out.ovl
  > diff ${expected_ovl} out/out.ovl

The new smart-hit-per-target. This should capture potential supplementary alignments of a query/target pair, but only on the same strand.
  $ rm -rf out && mkdir -p out
  > in_fasta=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.in.reads.plasmid_5kbp.ecoli_w_bc1087.fasta
  > expected_ovl=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.01.out.smart_hit_per_target.ovl
  > ${BIN_DIR}/pancake seqdb out/reads ${in_fasta}
  > ${BIN_DIR}/pancake seeddb out/reads.seqdb out/reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --num-threads 1 --skip-sym --write-rev --out-fmt ipa --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 out/reads out/reads 0 0 0 --out-fmt m4 --smart-hit-per-target --secondary-min-ovl-frac 0.05 | sort -k 1,1 > out/out.ovl
  > diff ${expected_ovl} out/out.ovl

#############################################
### Test the new 'smart' hits per target  ###
### On several datasets.                  ###
#############################################
Several plasmid reads from an E. Coli W dataset.
The new smart-hit-per-target. This should capture potential supplementary alignments of a query/target pair, but only on the same strand.
  $ rm -rf out && mkdir -p out
  > in_fasta=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.02.in.reads.plasmid_5kbp.ecoli_w_bc1092.fasta
  > expected_ovl=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.02.out.smart_hit_per_target.ovl
  > ${BIN_DIR}/pancake seqdb out/reads ${in_fasta}
  > ${BIN_DIR}/pancake seeddb out/reads.seqdb out/reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --num-threads 1 --skip-sym --write-rev --out-fmt ipa --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 out/reads out/reads 0 0 0 --out-fmt m4 --smart-hit-per-target --secondary-min-ovl-frac 0.05 | sort -k 1,1 > out/out.ovl
  > diff ${expected_ovl} out/out.ovl

Several plasmid reads from a K. Pneumoniae dataset.
The new smart-hit-per-target. This should capture potential supplementary alignments of a query/target pair, but only on the same strand.
  $ rm -rf out && mkdir -p out
  > in_fasta=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.03.in.reads.plasmid_2kbp.kpneumoniae_bc1074.fasta
  > expected_ovl=${TEST_DATA_DIR}/hifi-ovl/small-plasmids/test.03.out.smart_hit_per_target.ovl
  > ${BIN_DIR}/pancake seqdb out/reads ${in_fasta}
  > ${BIN_DIR}/pancake seeddb out/reads.seqdb out/reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --num-threads 1 --skip-sym --write-rev --out-fmt ipa --min-idt 98 --traceback --mask-hp --mask-repeats --trim --trim-window-size 30 --trim-match-frac 0.75 out/reads out/reads 0 0 0 --out-fmt m4 --smart-hit-per-target --secondary-min-ovl-frac 0.05 --min-map-len 200 --min-anchor-span 200 | sort -k 1,1 > out/out.ovl
  > diff ${expected_ovl} out/out.ovl

#############################################
### Sanity check that normal overlapping  ###
### works.                                ###
#############################################
Test a set of overlaps without supplementary same-hit alignments with the new "--smart-hit-per-target " to verify that it does not
filter excessively when supplementary hits do not exist.
Dovetail overlap 5prime test. Read "m64030_190330_071939/101844710/ccs" should have only 5prime overlaps with other reads.
This is a test from test_hifi_ovl.t.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile1-5prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --smart-hit-per-target | grep "^m64030_190330_071939/101844710/ccs"
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/62523930/ccs -7889 99.94 0 0 7895 11811 1 0 7894 11307 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61737670/ccs -6015 99.82 0 0 6031 11811 0 2564 8590 8590 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60884410/ccs -4942 99.68 0 0 4959 11811 0 4181 9139 9139 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61276980/ccs -3829 99.77 0 0 3845 11811 1 0 3838 10150 5
