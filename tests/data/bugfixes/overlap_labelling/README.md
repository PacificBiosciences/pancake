This small test dataset was extracted from a subsampled Drosophila dataset.
I noticed a strange case where the artificially generated symmetric overlap had a wrong label (3 prime instead of 5 prime):
```
pat_m64011_190605_003147/160434972/ccs pat_m64011_190605_003147/32900288/ccs -11736 99.660 0 0 11718 13995 1 0 11736 13966 3
pat_m64011_190605_003147/32900288/ccs pat_m64011_190605_003147/160434972/ccs -11736 99.660 0 0 11736 13966 1 0 11718 13995 5
```
The actual output should be something like:
```
pat_m64011_190605_003147/160434972/ccs pat_m64011_190605_003147/32900288/ccs -11736 99.660 0 0 11718 13995 1 0 11736 13966 5
pat_m64011_190605_003147/32900288/ccs pat_m64011_190605_003147/160434972/ccs -11736 99.660 0 0 11736 13966 1 0 11718 13995 5
```

This bug could cause problems for the chimera filter which parses this information.
For the string graph it wouldn't do much, because all overlaps are read and integrated into the graph.
