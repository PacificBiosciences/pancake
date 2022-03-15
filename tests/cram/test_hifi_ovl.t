Dovetail overlap 5prime test. Read "m64030_190330_071939/101844710/ccs" should have only 5prime overlaps with other reads.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile1-5prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/101844710/ccs"
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/62523930/ccs -7889 99.94 0 0 7895 11811 1 0 7894 11307 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61737670/ccs -6015 99.82 0 0 6031 11811 0 2564 8590 8590 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60884410/ccs -4942 99.68 0 0 4959 11811 0 4181 9139 9139 5
  m64030_190330_071939/101844710/ccs m64030_190330_071939/61276980/ccs -3829 99.77 0 0 3845 11811 1 0 3838 10150 5

Dovetail overlap 3prime test. Read "m64030_190330_071939/101909220/ccs" should have only 3prime overlaps with other reads.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile2-3prime.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/101909220/ccs"
  m64030_190330_071939/101909220/ccs m64030_190330_071939/165937510/ccs -7412 99.91 0 2111 9530 9530 1 2475 9898 9898 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/22611030/ccs -5388 99.15 0 4096 9530 9530 0 0 5435 9487 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/25626700/ccs -3280 99.91 0 6247 9530 9530 0 0 3286 10208 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/154862090/ccs -2988 99.57 0 6527 9530 9530 1 5043 8044 8044 3
  m64030_190330_071939/101909220/ccs m64030_190330_071939/157549180/ccs -1811 99.02 0 7696 9530 9530 1 8299 10128 10128 3

Dovetail overlap all fwd oriented.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -6944 99.27 0 0 6996 8602 0 656 7651 7651 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -6818 99.84 0 0 6829 8602 0 1987 8823 8823 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -6803 99.91 0 1793 8602 8602 0 0 6813 8913 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1783 99.78 0 0 1787 8602 0 7139 8928 8928 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -1324 99.85 0 7276 8602 8602 0 0 1328 9090 3

Dovetail overlap all rev oriented.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -9838 99.91 0 0 9847 10635 1 0 9851 10977 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -9406 99.87 0 1217 10635 10635 1 946 10371 10371 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -6167 99.85 0 4459 10635 10635 1 3541 9723 9723 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -4527 99.67 0 0 4547 10635 1 0 4542 9121 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -3645 99.75 0 6981 10635 10635 1 7165 10825 10825 3

Perfect matching overlap.
One sequence is a real CCS read ("m64030_190330_071939/101844710/ccs"), one is a fake read composed of 5000bp copy/pasted
from the real read, and the third sequence is the reverse complement of the fake read.
All should have perfect 100.00% identity matches.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

Test symmetry in the aligned identity when A->B and B->A.
TODO: This overlap is not symmetric because we could use trimming of the local alignments. It is context dependent how far
the alignment is extended in each direction.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile6-short-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --min-idt 0 reads reads 0 0 0
  m64030_190330_071939/39060280/ccs m64030_190330_071939/101975180/ccs -1929 94.36 0 5162 7255 9907 1 0 2047 10984 u
  m64030_190330_071939/101975180/ccs m64030_190330_071939/39060280/ccs -1929 94.38 0 0 2099 10984 1 5208 7255 9907 u

Contained overlap.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile7-contained.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/100008040/ccs m64030_190330_071939/117113640/ccs -9831 99.88 0 0 9847 9847 0 388 10231 11163 contained
  m64030_190330_071939/117113640/ccs m64030_190330_071939/100008040/ccs -9831 99.88 0 388 10231 11163 0 0 9847 9847 contains

Internal overlap.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile8-internal.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --min-idt 94.0 reads reads 0 0 0 | grep "^m64030_190330_071939/101451080/ccs"
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -5266 95.51 0 2168 7686 8022 0 5698 11212 11350 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -5232 95.48 0 2168 7648 8022 1 1238 6721 10056 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -4865 95.35 0 2043 7186 8022 1 0 5104 12064 u
  m64030_190330_071939/101451080/ccs m64030_190330_071939/160039920/ccs -1347 94.10 0 6263 7737 8022 1 7244 8678 8678 u
### Note2: these were the results with the SES-1986 algorithm. Internal overlaps don't have a clear boundary, so it's
### arbitrary how much we extend them beyond their "endpoint".
#  m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -3229 96.86 0 2140 5474 8022 1 1716 5063 12064 u
#  m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -3214 95.72 0 2116 5474 8022 0 5649 9017 11350 u
#  m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -3202 95.37 0 2116 5474 8022 1 3400 6770 10056 u
### Note: this differs from Raptor, Raptor is more sensitive here:
# m64030_190330_071939/101451080/ccs m64030_190330_071939/161024730/ccs -2633 96.7319 0 2256 7060 8022 1 126 4930 12064
# m64030_190330_071939/101451080/ccs m64030_190330_071939/164825240/ccs -2505 96.2606 0 2227 7415 8022 0 5753 10965 11350
# m64030_190330_071939/101451080/ccs m64030_190330_071939/164823320/ccs -2314 95.9715 0 2256 7444 8022 1 1428 6637 10056
# m64030_190330_071939/101451080/ccs m64030_190330_071939/165414370/ccs -1392 95.4594 0 4513 7390 8022 0 0 2797 10706
# m64030_190330_071939/101451080/ccs m64030_190330_071939/160039920/ccs -500 98.9011 0 6261 7357 8022 1 7586 8678 8678
### Also, to-do: larger '--aln-bw' actually can reduce sensitivity. Look into that.

Overlapping a single full pile.
(Awk is a hack to support cross-platform 'wc'. On some systems it outputs a prefix whitespace, and on some it does not.)
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile9-single-full-pile.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 | grep "^m64030_190330_071939/111673510/ccs"
  m64030_190330_071939/111673510/ccs m64030_190330_071939/93520790/ccs -10914 99.79 0 127 11064 11909 0 0 10953 10953 contains
  m64030_190330_071939/111673510/ccs m64030_190330_071939/150407580/ccs -10743 99.03 0 696 11544 11909 0 0 10858 10858 contains
  m64030_190330_071939/111673510/ccs m64030_190330_071939/69141390/ccs -10491 99.72 0 0 10520 11909 1 0 10540 11081 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/57868840/ccs -10023 99.79 0 1271 11315 11909 0 0 10058 10058 contains
  m64030_190330_071939/111673510/ccs m64030_190330_071939/93455110/ccs -9762 99.78 0 634 10418 11909 1 0 9799 9799 contains
  m64030_190330_071939/111673510/ccs m64030_190330_071939/15729630/ccs -9758 99.50 0 931 10738 11909 1 0 9830 9830 contains
  m64030_190330_071939/111673510/ccs m64030_190330_071939/54265740/ccs -8557 99.06 0 3271 11909 11909 0 0 8656 8870 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/155845340/ccs -8154 99.72 0 0 8177 11909 1 0 8193 9067 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/36505390/ccs -8137 99.60 0 2095 10265 11909 1 0 8181 8181 contains
  m64030_190330_071939/111673510/ccs m64030_190330_071939/68813360/ccs -8096 99.66 0 3785 11909 11909 1 1587 9714 9714 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/35586590/ccs -7996 98.71 0 0 8101 11909 1 0 8142 8354 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/115345600/ccs -7719 99.81 0 4175 11909 11909 0 0 7743 11204 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/7144350/ccs -7664 99.66 0 0 7690 11909 1 0 7703 11566 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/90049360/ccs -7475 99.76 0 0 7493 11909 1 0 7504 9876 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/153813000/ccs -7455 99.79 0 4438 11909 11909 1 4651 12132 12132 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/123994740/ccs -7454 99.68 0 0 7478 11909 1 0 7495 9495 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/93915010/ccs -7411 99.81 0 4484 11909 11909 0 0 7435 10541 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/50463970/ccs -7101 99.75 0 4790 11909 11909 0 0 7129 10843 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/8586840/ccs -6929 99.76 0 0 6946 11909 0 3166 10124 10124 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/24315480/ccs -6875 99.70 0 5013 11909 11909 1 2343 9246 9246 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/84674970/ccs -6841 99.30 0 0 6889 11909 0 2257 9160 9160 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/167381700/ccs -6454 99.77 0 0 6469 11909 0 2888 9367 9367 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/81921810/ccs -6396 99.80 0 5500 11909 11909 0 0 6418 9484 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/85591310/ccs -6395 99.52 0 0 6426 11909 1 0 6434 8478 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/158271420/ccs -6163 99.79 0 5733 11909 11909 1 3181 9366 9366 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/113706420/ccs -6091 99.72 0 0 6108 11909 0 3088 9206 9206 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/11010900/ccs -5970 99.04 0 0 6028 11909 1 0 6049 7599 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/58197860/ccs -5712 99.63 0 0 5733 11909 1 0 5741 8909 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/148507390/ccs -5586 99.57 0 6299 11909 11909 1 2485 8100 8100 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/53871380/ccs -4603 99.72 0 0 4616 11909 1 0 4624 7507 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/171051630/ccs -4446 99.71 0 0 4459 11909 1 0 4469 11305 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/60032210/ccs -4391 99.71 0 0 4404 11909 1 0 4414 8422 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/26936070/ccs -4301 99.06 0 7567 11909 11909 1 6433 10777 10777 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/119997480/ccs -4166 99.76 0 7733 11909 11909 1 6170 10354 10354 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/153682400/ccs -4119 98.83 0 7729 11909 11909 0 0 4168 12716 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/144967080/ccs -3444 99.37 0 0 3466 11909 0 5936 9411 9411 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/138150550/ccs -2968 99.57 0 8928 11909 11909 1 5759 8749 8749 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/21170930/ccs -2865 99.62 0 0 2876 11909 1 0 2885 9048 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/57871100/ccs -2749 99.93 0 9158 11909 11909 0 0 2753 9230 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/7602370/ccs -2327 99.91 0 9580 11909 11909 1 6457 8788 8788 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/94111340/ccs -2273 99.35 0 0 2288 11909 0 7881 10181 10181 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/50922520/ccs -2211 98.88 0 0 2236 11909 0 9000 11242 11242 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/122552600/ccs -2174 99.54 0 0 2184 11909 0 7535 9727 9727 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/163513680/ccs -2031 99.51 0 0 2041 11909 0 6436 8483 8483 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/44500860/ccs -1855 99.57 0 10043 11909 11909 1 8046 9909 9909 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/65733000/ccs -1757 99.55 0 0 1765 11909 1 0 1771 9936 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/40370270/ccs -1749 99.83 0 10157 11909 11909 1 7123 8878 8878 3
  m64030_190330_071939/111673510/ccs m64030_190330_071939/139331110/ccs -1666 99.52 0 0 1674 11909 1 0 1680 9562 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/92734730/ccs -1500 99.60 0 0 1506 11909 1 0 1510 7747 5
  m64030_190330_071939/111673510/ccs m64030_190330_071939/131140220/ccs -1403 99.43 0 0 1411 11909 0 9004 10419 10419 5

##################################################################
### Test removing symmetric arcs and writing reverse overlaps. ###
##################################################################
This input consists of a single dovetail overlap. By default, both arcs should be output.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Skip symmetric overlaps. The pair where Aid > Bid should be removed. (The ID is given in the DB, not the Aname and Bname.)
This input consists of a single dovetail overlap. By default, both arcs should be output.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Writing symmetric overlaps should provide the other arc of the pair of reads, even though only one of them was
actually computed. Order may change though.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym --write-rev reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
##################################################################

Test writing IDs instead of headers. Other than that, the same as the Perfect match overlap test.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --write-ids reads reads 0 0 0 --out-fmt m4
  000000000 000000001 -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  000000000 000000002 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  000000001 000000002 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  000000001 000000000 -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  000000002 000000001 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  000000002 000000000 -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

####################################################################
### Test allowing the output of self hits and symmetric overlaps ###
####################################################################
Normal run with default options. This should not output self hits and will not skip symmetric overlaps.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Self-hits allowed.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --allow-self-hits reads reads 0 0 0
  m64030_190330_071939/101844710/ccs m64030_190330_071939/101844710/ccs -11811 100.00 0 0 11811 11811 0 0 11811 11811 contained
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/60686570/ccs -10059 100.00 0 0 10059 10059 0 0 10059 10059 contained
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Self-hits NOT allowed, but different query and target DBs are used so Aid == Bid is allowed.
  $ ${BIN_DIR}/pancake seqdb reads1 ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads1.seqdb reads1
  > ${BIN_DIR}/pancake seqdb reads2 ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads2.seqdb reads2
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads1 reads2 0 0 0
  m64030_190330_071939/101844710/ccs m64030_190330_071939/101844710/ccs -11811 100.00 0 0 11811 11811 0 0 11811 11811 contained
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/60686570/ccs -10059 100.00 0 0 10059 10059 0 0 10059 10059 contained
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Skip symmetric overlaps.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym reads reads 0 0 0
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

Even though '--skip-sym' is specified, do not skip symmetric overlaps because different query and target DBs are provided.
Also, because the DBs are different, the quasi-"self" hits are allowed.
  $ ${BIN_DIR}/pancake seqdb reads1 ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads1.seqdb reads1
  > ${BIN_DIR}/pancake seqdb reads2 ${TEST_DATA_DIR}/hifi-ovl/reads.pile10-single-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads2.seqdb reads2
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --skip-sym reads1 reads2 0 0 0
  m64030_190330_071939/101844710/ccs m64030_190330_071939/101844710/ccs -11811 100.00 0 0 11811 11811 0 0 11811 11811 contained
  m64030_190330_071939/101844710/ccs m64030_190330_071939/60686570/ccs -9053 99.74 0 0 9077 11811 0 981 10059 10059 5
  m64030_190330_071939/60686570/ccs m64030_190330_071939/60686570/ccs -10059 100.00 0 0 10059 10059 0 0 10059 10059 contained
  m64030_190330_071939/60686570/ccs m64030_190330_071939/101844710/ccs -9053 99.74 0 981 10059 10059 0 0 9077 11811 3

##############
### Best N ###
##############
Bestn testing on a set of perfect matching overlaps. Value "0" should deactivate the filter.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn 0 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

Bestn testing on a set of perfect matching overlaps. Value < "0" should deactivate the filter.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn -1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5

Bestn testing on a set of perfect matching overlaps. Positive value should be applied normally.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn 1 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained

Bestn testing on a set of perfect matching overlaps. A large positive value, to test for potential edge cases.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile5-perfect-ovl.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 --bestn 10 reads reads 0 0 0 --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read/1/ccs -5000 100.00 0 0 5000 11811 0 10000 15000 15000 5
  m64030_190330_071939/101844710/ccs fake_read/1/ccs-inverted_0-15000 -5000 100.00 0 0 5000 11811 1 0 5000 15000 5
  fake_read/1/ccs fake_read/1/ccs-inverted_0-15000 -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs m64030_190330_071939/101844710/ccs -5000 100.00 0 10000 15000 15000 0 0 5000 11811 3
  fake_read/1/ccs-inverted_0-15000 fake_read/1/ccs -15000 100.00 0 0 15000 15000 1 0 15000 15000 contained
  fake_read/1/ccs-inverted_0-15000 m64030_190330_071939/101844710/ccs -5000 100.00 0 0 5000 15000 1 0 5000 11811 5
##############

###########################
### Traceback alignment ###
###########################
Traceback alignment.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -6968 99.27 0 0 6996 8602 0 656 7651 7651 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -6827 99.84 0 0 6829 8602 0 1987 8823 8823 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -6808 99.91 0 1793 8602 8602 0 0 6813 8913 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1786 99.78 0 0 1787 8602 0 7139 8928 8928 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -1326 99.85 0 7276 8602 8602 0 0 1328 9090 3

Traceback alignment.
Dovetail overlap all rev oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --out-fmt m4 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -9844 99.91 0 0 9847 10635 1 0 9851 10977 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -9415 99.87 0 1217 10635 10635 1 946 10371 10371 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -6174 99.85 0 4459 10635 10635 1 3541 9723 9723 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -4537 99.67 0 0 4547 10635 1 0 4542 9121 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -3652 99.75 0 6981 10635 10635 1 7165 10825 10825 3

Traceback and HPC turned on at the same time.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --use-hpc --min-map-len 0 --min-anchor-span 100 --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 99.90 0 0 5145 6298 0 496 5640 5640 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 99.98 0 0 5026 6298 0 1470 6497 6497 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3

Traceback and HPC turned on at the same time.
Dovetail overlap all rev oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile4-rev.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --use-hpc --min-map-len 0 --min-anchor-span 100 --out-fmt m4 | grep "^m64030_190330_071939/102172020/ccs"
  m64030_190330_071939/102172020/ccs m64030_190330_071939/43124430/ccs -7275 99.99 0 0 7276 7845 1 0 7276 8160 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/30016220/ccs -6944 99.99 0 900 7845 7845 1 705 7650 7650 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/28901470/ccs -4530 99.98 0 3314 7845 7845 1 2617 7148 7148 3
  m64030_190330_071939/102172020/ccs m64030_190330_071939/49610760/ccs -3379 99.97 0 0 3379 7845 1 0 3380 6838 5
  m64030_190330_071939/102172020/ccs m64030_190330_071939/52102340/ccs -2680 99.96 0 5164 7845 7845 1 5338 8019 8019 3

Writing the CIGAR string produced by the traceback. HPC is turned on at the same time just so that the CIGAR string is shorter.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Converts the input FASTA file into a SeqDB, computes the seeds using SeedDB, and then runs overlapping.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --use-hpc --min-map-len 0 --min-anchor-span 100 --write-cigar --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 99.90 0 0 5145 6298 0 496 5640 5640 5 221=1D649=1D296=1I169=1I167=1I3640=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 99.98 0 0 5026 6298 0 1470 6497 6497 5 5026=1D
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3 5008=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5 1288=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3 952=

Ignore indels when computing the alignment identity.
Dovetail overlap all fwd oriented, with the sensitive mode turned on.
Important: because of HPC compression, the overlap span can reduce below the default threshold of 1000bp. These
settings need to be modified below.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile3-fwd.fasta
  > ${BIN_DIR}/pancake seeddb --use-hpc reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --no-indels --use-hpc --min-map-len 0 --min-anchor-span 100 --write-cigar --out-fmt m4 | grep "^m64030_190330_071939/102303370/ccs"
  m64030_190330_071939/102303370/ccs m64030_190330_071939/108135170/ccs -5142 100.00 0 0 5145 6298 0 496 5640 5640 5 221=1D649=1D296=1I169=1I167=1I3640=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/109642940/ccs -5026 100.00 0 0 5026 6298 0 1470 6497 6497 5 5026=1D
  m64030_190330_071939/102303370/ccs m64030_190330_071939/106038710/ccs -5008 100.00 0 1290 6298 6298 0 0 5008 6585 3 5008=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/114034050/ccs -1288 100.00 0 0 1288 6298 0 5263 6551 6551 5 1288=
  m64030_190330_071939/102303370/ccs m64030_190330_071939/10946880/ccs -952 100.00 0 5346 6298 6298 0 0 952 6699 3 952=

Ignore SNPs when computing the alignment identity.
This is a synthetic test case.
Note: The mismatches introduced into the fake read are all near to homopolymers, so the SES aligner can easily
mistake them for an indel combination with interspersed matching bases.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile11-3_mismatches.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --min-map-len 0 --min-anchor-span 100 --write-cigar --no-snps
  m64030_190330_071939/101844710/ccs fake_read_with_3_snps/1/ccs -4997 100.00 0 0 5000 11811 0 10000 15000 15000 5 3000=1X999=1X898=1X100=
  fake_read_with_3_snps/1/ccs m64030_190330_071939/101844710/ccs -4997 100.00 0 10000 15000 15000 0 0 5000 11811 3 3000=1X999=1X898=1X100=

Ignore indels when computing the alignment identity.
This is a synthetic test case.
Note: The mismatches introduced into the fake read are all near to homopolymers, so the SES aligner can easily
mistake them for an indel combination with interspersed matching bases.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile11-3_mismatches.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --min-map-len 0 --min-anchor-span 100 --write-cigar --no-indels --out-fmt m4
  m64030_190330_071939/101844710/ccs fake_read_with_3_snps/1/ccs -4997 99.94 0 0 5000 11811 0 10000 15000 15000 5 3000=1X999=1X898=1X100=
  fake_read_with_3_snps/1/ccs m64030_190330_071939/101844710/ccs -4997 99.94 0 10000 15000 15000 0 0 5000 11811 3 3000=1X999=1X898=1X100=

#############################
### Test variant strings. ###
#############################
Overlap two reads where one should be reverse complemented.
This is a synthetic test case, where there is only one mismatches and the sequences are short.
The same alignment should be reported in fwd and rev directions, and it should be the same.
The variant strings should always be in the fwd strand of the sequence, so when the reads are swapped, the
actual variant sequences should in this case be the same (just swapped).
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  read1-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read4-rev read1-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *

An expansion of the previous test, but with more than 2 sequences.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile13-simple_errors_7_reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  read1-fwd read2-fwd -180 100.0000 0 0 180 180 0 0 180 180 c c u 180= * * *
  read1-fwd read1-rev -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read1-fwd read3-fwd -179 99.4444 0 0 180 180 0 0 179 179 c c u 54=1I125= G * *
  read1-fwd read4-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= T C *
  read1-fwd read3-rev -179 99.4444 0 0 180 180 1 0 179 179 c c u 54=1I125= G * *
  read1-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read2-fwd read1-fwd -180 100.0000 0 0 180 180 0 0 180 180 c c u 180= * * *
  read2-fwd read1-rev -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read2-fwd read3-fwd -179 99.4444 0 0 180 180 0 0 179 179 c c u 54=1I125= G * *
  read2-fwd read4-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= T C *
  read2-fwd read3-rev -179 99.4444 0 0 180 180 1 0 179 179 c c u 54=1I125= G * *
  read2-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read3-fwd read1-fwd -179 99.4444 0 0 179 179 0 0 180 180 c c u 54=1D125= * G *
  read3-fwd read2-fwd -179 99.4444 0 0 179 179 0 0 180 180 c c u 54=1D125= * G *
  read3-fwd read1-rev -179 99.4444 0 0 179 179 1 0 180 180 c c u 54=1D125= * C *
  read3-fwd read3-rev -179 100.0000 0 0 179 179 1 0 179 179 c c u 179= * * *
  read3-fwd read4-fwd -178 98.8889 0 0 179 179 0 0 180 180 c c u 54=1D97=1X27= T GC *
  read3-fwd read4-rev -178 98.8889 0 0 179 179 1 0 180 180 c c u 54=1D97=1X27= T GC *
  read4-fwd read4-rev -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read4-fwd read1-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= C T *
  read4-fwd read2-fwd -179 99.4444 0 0 180 180 0 0 180 180 c c u 152=1X27= C T *
  read4-fwd read1-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= C A *
  read4-fwd read3-fwd -178 98.8889 0 0 180 180 0 0 179 179 c c u 54=1I97=1X27= GC T *
  read4-fwd read3-rev -178 98.8889 0 0 180 180 1 0 179 179 c c u 54=1I97=1X27= GC A *
  read1-rev read1-fwd -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read1-rev read2-fwd -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read1-rev read3-fwd -179 99.4444 0 0 180 180 1 0 179 179 c c u 125=1I54= C * *
  read1-rev read4-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= A C *
  read1-rev read3-rev -179 99.4444 0 0 180 180 0 0 179 179 c c u 125=1I54= C * *
  read1-rev read4-rev -179 99.4444 0 0 180 180 0 0 180 180 c c u 27=1X152= A G *
  read3-rev read1-fwd -179 99.4444 0 0 179 179 1 0 180 180 c c u 125=1D54= * G *
  read3-rev read2-fwd -179 99.4444 0 0 179 179 1 0 180 180 c c u 125=1D54= * G *
  read3-rev read3-fwd -179 100.0000 0 0 179 179 1 0 179 179 c c u 179= * * *
  read3-rev read1-rev -179 99.4444 0 0 179 179 0 0 180 180 c c u 125=1D54= * C *
  read3-rev read4-fwd -178 98.8889 0 0 179 179 1 0 180 180 c c u 27=1X97=1D54= A GC *
  read3-rev read4-rev -178 98.8889 0 0 179 179 0 0 180 180 c c u 27=1X97=1D54= A GC *
  read4-rev read4-fwd -180 100.0000 0 0 180 180 1 0 180 180 c c u 180= * * *
  read4-rev read1-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *
  read4-rev read2-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *
  read4-rev read1-rev -179 99.4444 0 0 180 180 0 0 180 180 c c u 27=1X152= G A *
  read4-rev read3-fwd -178 98.8889 0 0 180 180 1 0 179 179 c c u 27=1X97=1I54= GC T *
  read4-rev read3-rev -178 98.8889 0 0 180 180 0 0 179 179 c c u 27=1X97=1I54= GC A *

#########################################
### Test writing in different formats ###
#########################################
M4 output format - with traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt m4
  read1-fwd read4-rev -179 99.44 0 0 180 180 1 0 180 180 contained 152=1X27=
  read4-rev read1-fwd -179 99.44 0 0 180 180 1 0 180 180 contained 27=1X152=

IPAOvl output format - with traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  read1-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  read4-rev read1-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *

SAM output format - with traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt sam
  @HD	VN:1.5
  @SQ	SN:read1-fwd	LN:180
  @SQ	SN:read4-rev	LN:180
  read1-fwd	16	read4-rev	1	60	27=1X152=	*	0	0	ATTATTGGCTCACTCTTACTCAGCCTTATTGTGGCAATGGTCATGAAGCGCAATAAAACCGTATAACATCTCTCTGTGCGCAGTACTTCCTGTATTATTGTGGTGGCGGTCGATATGCGCACTGGCAAAAAAACGGGCTTGAATATCGTTGAAACCCTTTAACAAAGCACAGGAGCGTGC	*	AT:Z:c	BT:Z:c
  read4-rev	16	read1-fwd	1	60	152=1X27=	*	0	0	GCACGCTCCTGTGCTTTGTTAAAGGGTTTCAACGATATTCAAGCCCGTTTTTTTGCCAGTGCGCATATCGACCGCCACCACAATAATACAGGAAGTACTGCGCACAGAGAGATGTTATACGGTTTTATTGCGCTTCATGACCATTGCCACAACAAGGCTGAGTAAGAGTGAGCCAATAAT	*	AT:Z:c	BT:Z:c

PAF output format - with traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt paf
  read1-fwd	180	0	180	-	read4-rev	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	cg:Z:152=1X27=	VQ:Z:T	VT:Z:G
  read4-rev	180	0	180	-	read1-fwd	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	cg:Z:27=1X152=	VQ:Z:G	VT:Z:T



M4 output format - no traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt m4
  read1-fwd read4-rev -179 99.44 0 0 180 180 1 0 180 180 contained
  read4-rev read1-fwd -179 99.44 0 0 180 180 1 0 180 180 contained

IPAOvl output format - no traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  read1-fwd read4-rev -179 99.4444 0 0 180 180 1 0 180 180 c c u * * * *
  read4-rev read1-fwd -179 99.4444 0 0 180 180 1 0 180 180 c c u * * * *

SAM output format - no traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt sam
  @HD	VN:1.5
  @SQ	SN:read1-fwd	LN:180
  @SQ	SN:read4-rev	LN:180
  read1-fwd	16	read4-rev	1	60	*	*	0	0	ATTATTGGCTCACTCTTACTCAGCCTTATTGTGGCAATGGTCATGAAGCGCAATAAAACCGTATAACATCTCTCTGTGCGCAGTACTTCCTGTATTATTGTGGTGGCGGTCGATATGCGCACTGGCAAAAAAACGGGCTTGAATATCGTTGAAACCCTTTAACAAAGCACAGGAGCGTGC	*	AT:Z:c	BT:Z:c
  read4-rev	16	read1-fwd	1	60	*	*	0	0	GCACGCTCCTGTGCTTTGTTAAAGGGTTTCAACGATATTCAAGCCCGTTTTTTTGCCAGTGCGCATATCGACCGCCACCACAATAATACAGGAAGTACTGCGCACAGAGAGATGTTATACGGTTTTATTGCGCTTCATGACCATTGCCACAACAAGGCTGAGTAAGAGTGAGCCAATAAT	*	AT:Z:c	BT:Z:c

PAF output format - no traceback.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt paf
  read1-fwd	180	0	180	-	read4-rev	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	VQ:Z:*	VT:Z:*
  read4-rev	180	0	180	-	read1-fwd	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	VQ:Z:*	VT:Z:*


M4 output format - with traceback. Writing IDs.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --write-ids --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt m4
  000000000 000000001 -179 99.44 0 0 180 180 1 0 180 180 contained 152=1X27=
  000000001 000000000 -179 99.44 0 0 180 180 1 0 180 180 contained 27=1X152=

IPAOvl output format - with traceback. Writing IDs.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --write-ids --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt ipa
  000000000 000000001 -179 99.4444 0 0 180 180 1 0 180 180 c c u 152=1X27= T G *
  000000001 000000000 -179 99.4444 0 0 180 180 1 0 180 180 c c u 27=1X152= G T *

SAM output format - with traceback. Writing IDs.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --write-ids --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt sam
  @HD	VN:1.5
  @SQ	SN:000000000	LN:180
  @SQ	SN:000000001	LN:180
  000000000	16	000000001	1	60	27=1X152=	*	0	0	ATTATTGGCTCACTCTTACTCAGCCTTATTGTGGCAATGGTCATGAAGCGCAATAAAACCGTATAACATCTCTCTGTGCGCAGTACTTCCTGTATTATTGTGGTGGCGGTCGATATGCGCACTGGCAAAAAAACGGGCTTGAATATCGTTGAAACCCTTTAACAAAGCACAGGAGCGTGC	*	AT:Z:c	BT:Z:c
  000000001	16	000000000	1	60	152=1X27=	*	0	0	GCACGCTCCTGTGCTTTGTTAAAGGGTTTCAACGATATTCAAGCCCGTTTTTTTGCCAGTGCGCATATCGACCGCCACCACAATAATACAGGAAGTACTGCGCACAGAGAGATGTTATACGGTTTTATTGCGCTTCATGACCATTGCCACAACAAGGCTGAGTAAGAGTGAGCCAATAAT	*	AT:Z:c	BT:Z:c

PAF output format - with traceback. Writing IDs.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --write-ids --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt paf
  000000000	180	0	180	-	000000001	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	cg:Z:152=1X27=	VQ:Z:T	VT:Z:G
  000000001	180	0	180	-	000000000	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	cg:Z:27=1X152=	VQ:Z:G	VT:Z:T

PAF output format, generating symetric overlaps - with traceback. Writing IDs.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/hifi-ovl/reads.pile12-simple_errors_2-reads.fasta
  > ${BIN_DIR}/pancake seeddb -k 15 -w 10 -s 0 reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --write-ids --traceback --write-cigar --min-map-len 0 --min-anchor-span 20 --aln-bw 0.10 --aln-diff-rate 0.10 --out-fmt paf --skip-sym --write-rev
  000000001	180	0	180	-	000000000	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	cg:Z:27=1X152=	VQ:Z:G	VT:Z:T
  000000000	180	0	180	-	000000001	180	0	180	180	180	60	tp:A:P	NM:i:1	IT:f:99.4444	SC:i:-179	AT:Z:c	BT:Z:c	cg:Z:152=1X27=	VQ:Z:T	VT:Z:G

##########################################
### Testing masking of variant strings ###
##########################################
Write the variant strings and mask homopolymers and simple repeats. Align only one arc.
I manually verified the variants in this alignment using IGV.
This is the list of variants broken down by their positions (note that Bvars is reverse complemented).
OP      END     A   B       Valid   Description
385=    385
1I      1622    t           +       HP
1237=   1623
1I      1624    t           +       HP
1=      1625
1I      1626    T           +
930=    2556
1I      2557    c           +       HP
75=     2632
1D      2632        t       +       HP
57=     2689
1X      2690    T   G       +       This SNP is in a HP of 3 T bases, but it is not an indel so not masked.
8=      2698
1D      2698        g       +       HP
61=     2759
1I      2760    G           +
780=    3540
1I      3541    g           +       HP
680=    4221
1D      4221        c       +       HP
165=    4386
1D      4386        t       +       HP
545=    4931
1I      4932    c           +       HP
442=    5374
1I      5375    c           +       HP
730=    6105
1I      6106    t           +       HP
589=    6695
1D      6695        G       +       This G is placed in between two C bases in 000000551. Technically not a homopolymer error?
370=    7065
1D      7065        g       +       HP
119=    7184
1I      7185    G           +       This G is placed in between two C bases in 000000571. Technically not a homopolymer error?
463=    7648
1D      7648        T       +
639=    8287
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads1.weird_masking_issue.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --mask-hp --mask-repeats --skip-sym --write-rev --out-fmt ipa
  000000571 000000551 -8276 99.9275 0 0 8287 10491 1 0 8284 8714 5 5 u 385=1I1237=1I1=1I930=1I75=1D57=1X8=1D61=1I780=1I680=1D165=1D545=1I442=1I730=1I589=1D370=1D119=1I463=1D639= ttTcTGgcctG AcCagcCa *
  000000551 000000571 -8276 99.9275 0 0 8284 8714 1 0 8287 10491 5 5 u 639=1I463=1D121=1I368=1I590=1D730=1D443=1D544=1I165=1I680=1D779=1D64=1I5=1X58=1I76=1D928=1D2=1D1238=1D383= AcCagcCa ttTcTGgcctG *

Exactly the same test as one above, but the HP SNP is masked so the identity is slightly higher.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads1.weird_masking_issue.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --mask-hp --mask-repeats --mask-hp-snps --skip-sym --write-rev --out-fmt ipa
  000000571 000000551 -8276 99.9396 0 0 8287 10491 1 0 8284 8714 5 5 u 385=1I1237=1I1=1I930=1I75=1D57=1X8=1D61=1I780=1I680=1D165=1D545=1I442=1I730=1I589=1D370=1D119=1I463=1D639= ttTctGgcctG AcCagcca *
  000000551 000000571 -8276 99.9396 0 0 8284 8714 1 0 8287 10491 5 5 u 639=1I463=1D121=1I368=1I590=1D730=1D443=1D544=1I165=1I680=1D779=1D64=1I5=1X58=1I76=1D928=1D2=1D1238=1D383= AcCagcca ttTctGgcctG *

Also exactly the same as 2 tests above, but with masking of arbitrary HP indels.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads1.weird_masking_issue.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --mask-hp --mask-hp-arbitrary --mask-repeats --skip-sym --write-rev --out-fmt ipa
  000000571 000000551 -8276 99.9517 0 0 8287 10491 1 0 8284 8714 5 5 u 385=1I1237=1I1=1I930=1I75=1D57=1X8=1D61=1I780=1I680=1D165=1D545=1I442=1I730=1I589=1D370=1D119=1I463=1D639= ttTcTGgcctg AccagcCa *
  000000551 000000571 -8276 99.9517 0 0 8284 8714 1 0 8287 10491 5 5 u 639=1I463=1D121=1I368=1I590=1D730=1D443=1D544=1I165=1I680=1D779=1D64=1I5=1X58=1I76=1D928=1D2=1D1238=1D383= AccagcCa ttTcTGgcctg *


Test identity recomputation after masking.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads3.masking_and_identity.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --mask-hp --mask-repeats --skip-sym --write-rev --out-fmt ipa --min-idt 99.0
  000258100 000000002 -9511 100.0000 0 0 9514 10599 1 0 9512 13810 5 5 u 667=1D4901=1I3385=1I257=1I301= tga a *
  000000002 000258100 -9511 100.0000 0 0 9512 13810 1 0 9514 10599 5 5 u 303=1D256=1D3388=1D4902=1I662= a tga *

Test masking of indels at the very beginning, when the front of the sequence is a HP.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads4.mask_indels_at_front.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --mask-hp --mask-repeats --skip-sym --write-rev --out-fmt ipa
  001235617 000125429 -6445 99.9380 0 0 6461 21031 0 12793 19250 19250 5 3 u 1=1D758=1D103=1I158=1X87=1D304=1I73=1D32=1D323=1I84=1I123=1I98=1I175=1D204=1I12=1D5=1I64=1I61=1D54=1X54=1I1150=1D807=1I1416=1I98=1I94=1I53=1D54= cAggtcccgaTaTaat TtTtataaaAtt *
  000125429 001235617 -6445 99.9380 0 12793 19250 19250 0 0 6461 21031 3 5 u 1=1I758=1I103=1D158=1X87=1I304=1D73=1I32=1I323=1D84=1D123=1D98=1D175=1I204=1D12=1I5=1D64=1D61=1I54=1X54=1D1150=1I807=1D1416=1D98=1D94=1D53=1I54= TtTtataaaAtt cAggtcccgaTaTaat *

Test masking of indels at the very beginning, when the front of the sequence is a HP. Same as before, but mask arbitrary indels in HPs as well.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads4.mask_indels_at_front.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --write-cigar --mask-hp --mask-hp-arbitrary --mask-repeats --skip-sym --write-rev --out-fmt ipa
  001235617 000125429 -6445 99.9690 0 0 6461 21031 0 12793 19250 19250 5 3 u 1=1D758=1D103=1I158=1X87=1D304=1I73=1D32=1D323=1I84=1I123=1I98=1I175=1D204=1I12=1D5=1I64=1I61=1D54=1X54=1I1150=1D807=1I1416=1I98=1I94=1I53=1D54= cAggtcccgaTataat ttTtataaaAtt *
  000125429 001235617 -6445 99.9690 0 12793 19250 19250 0 0 6461 21031 3 5 u 1=1I758=1I103=1D158=1X87=1I304=1D73=1I32=1I323=1D84=1D123=1D98=1D175=1I204=1D12=1I5=1D64=1D61=1I54=1X54=1D1150=1I807=1D1416=1D98=1D94=1D53=1I54= ttTtataaaAtt cAggtcccgaTataat *

Same test at the first one in this section, but the output is to be written without the CIGAR and variant strings.
  $ ${BIN_DIR}/pancake seqdb reads ${TEST_DATA_DIR}/varstr/reads1.weird_masking_issue.fasta
  > ${BIN_DIR}/pancake seeddb reads.seqdb reads
  > ${BIN_DIR}/pancake ovl-hifi --num-threads 1 reads reads 0 0 0 --traceback --mask-hp --mask-repeats --skip-sym --write-rev --out-fmt ipa
  000000571 000000551 -8276 99.9275 0 0 8287 10491 1 0 8284 8714 5 5 u * * * *
  000000551 000000571 -8276 99.9275 0 0 8284 8714 1 0 8287 10491 5 5 u * * * *
