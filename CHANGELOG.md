# Changelog - Pancake

## Active version in development
### Changes
 - Add `dynamicBandwidth` boolean to `AlignmentParameters` that allows iterative KSW2 bandwidth exploration

## v1.3.0
### Version
- `pancake` - commit 1c1b037270bdd39b96aafc4e537345878cf4c923 (April 02, 2021), `pancake 1.3.0 (commit SL-release-10.0.0-557-g1c1b037)`

### Changes
- Arbitrary sequence IDs - sequences can now have an arbitrary numeric ID for mapping. Previously, the IDs of sequences needed to be consecutive and starting from 0.
- CollectHits now reverses the query coordinate instead of the target coordinate if the seeds are on the opposite strands. This allows for simplification of a lot of the code, and perhaps some efficiency improvement.
- Minor refactoring of the SeedIndex.
- GPU CIGAR conversion for the MapperBatchGPU. This is no longer performed on the CPU, and this improves the efficiency of this mapper.
- Removed the `maxHPCLen` parameter from minimizer computation. There is no longer an upper limit on the length of the homopolymers for compression. The `--max-hpc-len` parameter is now removed from `pancake seeddb`.
- Bugfixes for edge cases in minimizer computation: first window had a bug when multiple equal keys are there; spaced seeds could have omitted minimizers in the first window.
- Bugfix for negative seed hit spans (when the span of a seed ran over the range of the signed 8-bit ints).

## v1.2.0
### Version
- `pancake` - commit 2f05f3e16219291485de130f7ca68dc40ae7c5c5 (Mar 11, 2021), `pancake 1.2.0 (commit SL-release-10.0.0-417-g2f05f3e)`

### Changes
- The `AlignmentResult` now stores the number of alignment differences, so that the client doesn't have to parse the CIGAR string again. Updated the aligners to produce the diff counts.
- `FastaSequenceCachedStore` implementation. This is a store for `FastaSequenceCached` objects, which are like a "view" into the FASTA sequences (they do not own the data, only point to it). This is now used in the `SeqDBReaderCachedBlock`.
- Generic interface to `MapperHiFi`. It's now possible to run a single function which takes a set of `std::string` targets and queries, and maps/aligns them.
- Added the SeqDB dump tool (`seqdb-dump`) to dump the entire DB or just one block.
- Added the SeqDB info tool (`seqdb-info`) to display information about the SeqDB (length, N50, etc.).
- Minor refactoring, versioning code.
- (non-IPA related) KSW2 third-party library now removed from Pancake, and used from Pbcopper.
- (non-IPA related) Refactoring the MapperCLR: parametrized the flankExtensionFactor setting ; and always computing alignment regions during mapping, instead of only when performing alignment.
- (non-IPA related) `AlignerBatchCPU` implementation - an aligner which first takes a batch of sequence pairs, then performs alignment of all submitted pairs in parallel, on CPU.
- (non-IPA related) `MapperBatchCPU` implementation - a batch mapper which first takes a batch of sequence pairs, then performs mapping of sequences in parallel on the CPU. Following that, the mappings are aligned in batch using the AlignerBatchCPU in parallel on the CPU.
- (non-IPA related) `AlignerBatchGPU` - analogous to the AlignerBatchGPU, but uses the Cudaaligner (GenomeWorks) to perform alignment of the sequence pairs on the GPU.
- (non-IPA related) `MapperBatchGPU` - analogous to the `MapperBatchCPU`, but this uses the GPU aligner to align sequences. Mapping is still performed on the CPU, same as with the `MapperBatchCPU`.
- (non-IPA related) Abstract class `MapperBatchBase` added, and `MapperBatchCPU` and `MapperBatchGPU` now both inherit the interface.
- (non-IPA related) Optional dependency injection of the FireAndForget object into the `AlignerBatchCPU`, `MapperBatchCPU` and `MapperBatchGPU` to better handle kernel overload.
- All CLR mappers (MapperCLR, MapperBatchCPU, MapperBatch GPU) can now optionally skip/mock perfect alignments, and skip symmetric overlaps. (Mocking means that if a self-hit is detected based on sequence IDs, a perfect alignment result will be generated without actually computing the alignment.)
- Fixed the bestNSecondary in MapperCLR::Map_ where it did not properly respect the value < 0. (Instead of reporting all alignments, only some were reported.)
- The batch mappers (MapperBatchCPU, MapperBatchGPU) can accept custom mapping parameters for every chunk in a batch.
- Converted the `unique_ptr` to references in the overlap writers.
- Fixed the score computation in MapperCLR.

## v1.1.0 - SL - Release 10.1.0
### Version
- `pancake` - commit 29844f96d5cf58874a04ffd6a9895abe00a0750b (origin/release/prep) (Dec 16, 2020), `pancake 1.1.0 (commit SL-release-10.0.0-186-g29844f9)`

### Changes
- LIS chaining of seed hits. Improvements in overlap sensitivity by refining seed hits before alignment. The seed hits are now refined using the Longest Increasing Subsequence (LIS) algorithm.
- Scalability improvement in Pancake - reduced system time by removing `std::clock()` usage. It turns out that this function is blocking.
- Filtering duplicate overlaps in Pancake, which can occur when there is a larger gap in between seed hits.
- Pbcopper subproject now points to the GitHub mirror.
- GitHub mirror of Pancake now exists: https://github.com/PacificBiosciences/pancake
- Handling edge cases in alignment explicitly (e.g. when either query or target is of length `0`). Some aligners behaved non-deterministically in these cases.
- Resolved warnings.
- Minor refactoring.
- Encoded seed span. Seed span is now encoded in the Seed structure. The same as in Minimap2. Span of a seed hit in both query and target is now encoded in the SeedHit structure. Seed indexing and hit collection now works correctly in the homopolymer-compressed mode. It also works with the spaced seeds feature. The new maximum kmer size is equal to `28` now. This is because the span is encoded within the seed data structure, so some bits had to be specialized. A lot of tests had to be updated because previously they used the default of `k = 30`.
- (non-IPA related) Added MapperCLR - implemented a CLR mapper, currently available via API.  and several more followups to resolve warnings and improve API

## v0.2.0 - First release - SL - Release 10.0.0 (SL-release-10.0.0)
### Version
- `pancake` - commit dd98868876ff0171be3e191973834e1ac5229123 (HEAD, tag: SL-release-10.0.0, origin/release/10.0.0) (Sep 23, 2020) `pancake 0.2.0 (commit SL-release-10.0.0)`
