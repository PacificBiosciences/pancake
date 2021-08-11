/*
 * Minimizers.cpp
 *
 * This source file was obtained and adjusted from the Raptor
 * graph-based mapping tool codebase.
 *
 *  Created on: Oct 04, 2017
 *      Author: Ivan Sovic
 */

#include <pacbio/pancake/Lookups.h>
#include <pacbio/pancake/Minimizers.h>
#include <pacbio/pancake/Seed.h>
#include <pacbio/util/CommonTypes.h>
#include <deque>
#include <iostream>

namespace PacBio {
namespace Pancake {
namespace SeedDB {

using minkey_t = uint64_t;

const int32_t MAX_SPACING_IN_SEED = 32;
const int32_t MAX_WINDOW_BUFFER_SIZE = 512;

// #define DEBUG_GENERATE_MINIMIZERS

// This is a temporary class, just used for implementing the minimizer window.
class SpacedBuffer
{
public:
    SpacedBuffer()
    {
        for (size_t i = 0; i < seqBuffer.size(); ++i) {
            seqBuffer[i] = 0x0;
            seqBufferRC[i] = 0x0;
            numBasesIn[i] = 0;
        }
    }

    ~SpacedBuffer() = default;

    // Holds the current 2-bit seed representation.
    std::array<uint64_t, MAX_SPACING_IN_SEED> seqBuffer;
    // Holds the reverse complement 2-bit seed at the same position.
    std::array<uint64_t, MAX_SPACING_IN_SEED> seqBufferRC;
    // Number of bases added to the buffer.
    std::array<int32_t, MAX_SPACING_IN_SEED> numBasesIn;

    // Define the new circular buffer for the window.
    Seed winBuff[MAX_WINDOW_BUFFER_SIZE];

    int32_t winBuffPos = 0;
    int32_t winBuffMinPos = -1;
    int32_t kmerSpan = 0;
    std::deque<int32_t> hpEvents;

    void Clear(int32_t winSize)
    {
        // buffer = 0;
        // buffer_rc = 0;
        // numBasesIn = 0;
        winBuffPos = 0;
        winBuffMinPos = -1;
        kmerSpan = 0;
        hpEvents.clear();
        for (int32_t i = 0; i < winSize; ++i) {
            winBuff[i] = Seed();
        }
        for (size_t i = 0; i < seqBuffer.size(); ++i) {
            seqBuffer[i] = 0x0;
            seqBufferRC[i] = 0x0;
            numBasesIn[i] = 0;
        }
    }
};

int GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& minimizers, const uint8_t* seq,
                       const int32_t seqLen, const int32_t seqOffset, const int32_t seqId,
                       const int32_t kmerSize, const int32_t winSize, const int32_t spacing,
                       const bool useReverseComplement, const bool useHPC)
{
    // clang-format off

    const int32_t MAX_SEED_SPAN = 255;

    // Sanity check that the seq is not NULL;
    if (!seq) {
        throw std::runtime_error("Cannot generate minimizers. The sequence is NULL.");
    }
    // Sanity check for the size of kmer. It can't
    // be too large, or it won't fit the uint64_t buffer.
    if (kmerSize <= 0 || kmerSize > 28) {
        throw std::runtime_error(
            "Cannot generate minimizers. The kmerSize is out of bounds, should be in range [0, "
            "28]. kmerSize = " +
            std::to_string(kmerSize));
    }
    if (winSize > MAX_WINDOW_BUFFER_SIZE) {
        throw std::runtime_error(
            "Cannot generate minimizers. The winSize is out of bounds, should be in range [0, " +
            std::to_string(MAX_WINDOW_BUFFER_SIZE) + "]. winSize = " + std::to_string(winSize));
    }
    if (spacing > MAX_SPACING_IN_SEED) {
        throw std::runtime_error(
            "Cannot generate minimizers. The spacing is out of bounds, should be in range [0, " +
            std::to_string(MAX_SPACING_IN_SEED) + "]. spacing = " + std::to_string(spacing));
    }
    // Not technically an error if the seqLen is 0,
    // but there's nothing to do, so return.
    if (seqLen == 0) {
        return 0;
    }

    #ifdef DEBUG_GENERATE_MINIMIZERS
        std::cerr << "[GenerateMinimizers] seqId = " << seqId << ", seqLen = " << seqLen << "\n";
    #endif

    // Mask the number of required bits for the kmer.
    const minkey_t mask = (kmerSize < 32) ? ((((uint64_t)1) << (2 * kmerSize)) - 1) : 0xFFFFFFFFFFFFFFFF;

    SpacedBuffer bufferWD;
    const int32_t minSeedSpanSize = kmerSize * (spacing + 1) - spacing; // Equivalent to: `k + (k - 1) * s`

    auto IsAlternativeMinimizer = [](const Seed& s1, const Seed& s2) -> bool {
        // Another seed is a potential alternative minimizer only if it has
        // the same key and span but different position.
        return (s1.key == s2.key && s1.span == s2.span && (s1.pos != s2.pos || s1.seqID != s2.seqID || s1.seqRev != s2.seqRev));
    };

    auto IsWindowBufferElementValid = [](const Seed* winBuff, const int32_t elementId) -> bool {
        if (elementId < 0 || elementId >= MAX_WINDOW_BUFFER_SIZE) {
            return false;
        }
        return winBuff[elementId].Valid();
    };

    auto IsWindowFull = [](const int32_t numBasesInForSpace, const int32_t _winSize, const int32_t _kmerSize, const int32_t _spacing, const int32_t minSeedSpan, const int32_t offset) -> bool {
        const int32_t currentSpan = numBasesInForSpace + (numBasesInForSpace - 1) * _spacing;
        const bool ret = (numBasesInForSpace >= _kmerSize && currentSpan >= (_winSize + minSeedSpan + offset));
        #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "    [IsWindowFull] result = " << ret << "; numBasesInForSpace = " << numBasesInForSpace << ", currentSpan = " << currentSpan << ", _winSize = " << _winSize
                << ", _kmerSize = " << _kmerSize << ", _spacing = " << _spacing << ", minSeedSpan = " << minSeedSpan << ", offset = " << offset << "\n";
        #endif
        return ret;
    };

    auto IsWindowFullFirstTime = [](const int32_t numBasesInForSpace, const int32_t _winSize, const int32_t _kmerSize, const int32_t _spacing, const int32_t minSeedSpan) -> bool {
        const int32_t currentSpan = numBasesInForSpace + (numBasesInForSpace - 1) * _spacing;
        const bool ret = (numBasesInForSpace >= _kmerSize && currentSpan == (_winSize + minSeedSpan - 1));
        #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "    [IsWindowFullFirstTime] result = " << ret << "; numBasesInForSpace = " << numBasesInForSpace << ", currentSpan = " << currentSpan << ", _winSize = " << _winSize
                << ", _kmerSize = " << _kmerSize << ", _spacing = " << _spacing << ", minSeedSpan = " << minSeedSpan << "\n";
        #endif
        return ret;
    };

    for (int32_t pos = 0, space = 0; pos < seqLen; ++pos, ++space) {
        const int8_t b = seq[pos];
        Seed newSeed;

        // The new seed is automatically not valid, but let's be future-proof.
        newSeed.SetInvalid();

        if (space > spacing) {
            space = 0;
        }

        #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "[pos = " << pos << ", space = " << space << "; k = " << kmerSize << ", w = " << winSize << ", s = " << spacing << "]"
                            << " State: b = '" << b << "', winBuffPos = " << bufferWD.winBuffPos
                            << ", numBasesIn[space] = " << bufferWD.numBasesIn[space]
                            << ", winBuffMinPos = " << bufferWD.winBuffMinPos
                            << ", minimizer = [" << (bufferWD.winBuffMinPos < 0 ? "(no_minimizer)" : bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose()) << "]"
                            << "\n";
        #endif

        if (IsNucleotide[b]) {
            // If enabled, skip same bases and compute the span.
            if (useHPC) {
                // In this case, find the stretch of homopolymers (at least 1 base), and
                // add that stretch to "seed_span" (the seed_span is iteratively increased every
                // loop iteration).
                int32_t hpLen = 1;
                for (hpLen = 1; (pos + hpLen) < seqLen; ++hpLen) {
                    if (seq[pos + hpLen] != b) {
                        break;
                    }
                }
                bufferWD.hpEvents.push_back(hpLen);
                pos += hpLen - 1;
                bufferWD.kmerSpan += hpLen;
                // We need to keep track of the length of each HP event. They contribute as a single
                // base, but in reality, the span is larger. This deque keeps track of all hp events
                // in the current kmer.
                while (bufferWD.hpEvents.size() > static_cast<size_t>(kmerSize)) {
                    bufferWD.kmerSpan -= bufferWD.hpEvents.front();
                    bufferWD.hpEvents.pop_front();
                }
            } else {
                bufferWD.kmerSpan = std::min((bufferWD.numBasesIn[space] + 1) * (1 + spacing) - spacing, minSeedSpanSize);
            }

            // Add the base to the buffer.
            bufferWD.seqBuffer[space] = ((bufferWD.seqBuffer[space] << 2) | (((static_cast<uint64_t>(BaseToTwobit[b]))))) & mask;
            bufferWD.seqBufferRC[space] = ((bufferWD.seqBufferRC[space] >> 2) | (((static_cast<uint64_t>(BaseToTwobitComplement[b]))) << (kmerSize * 2 - 2))) & mask;

            // Determine the orientation of the key.
            minkey_t key = bufferWD.seqBuffer[space];
            bool isRev = false;
            if (useReverseComplement && bufferWD.seqBufferRC[space] < key) {
                key = bufferWD.seqBufferRC[space];
                isRev = true;
            }

            ++bufferWD.numBasesIn[space];

            if (bufferWD.numBasesIn[space] >= kmerSize) {
                // The 'pos' is the current position which is inclusive. We need to add a +1 to make it non-inclusive, so that the start position is calculated properly.
                // The kmerStart has to be computed here because we don't know the current kmerSpan yet (HP).
                const int32_t kmerStart = (pos + 1) - (bufferWD.kmerSpan);
                newSeed = Seed(InvertibleHash(key, mask), bufferWD.kmerSpan, seqId, kmerStart + seqOffset, isRev);
                // This is not necessary because the constructor will mark invalid
                // spans automatically, but let's make sure here and be future-proof.
                if (bufferWD.kmerSpan > MAX_SEED_SPAN) {
                    newSeed.SetInvalid();
                }

                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "    -> Constructed: " << newSeed << "\n";
                #endif
            }

            // Skip symmetric kmers.
            if (bufferWD.seqBuffer[space] == bufferWD.seqBufferRC[space]) {
                newSeed.SetInvalid();
            }

        } else {
            // If we encountered a non-nucleotide base, this is enough to reset the minimizer coding.
            // First, write out the current minimizer so it doesn't get lost in the void.
            // We only need to write one and not loop through all of them, because when we encounter a seed
            // with the same key then we write out the previous minimizer out. This means that we only need
            // to write out the current minimizer.

            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    - Resetting buffers.\n";
            #endif

            if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos) && IsWindowFull(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize, 0)) {
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());

                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "    - Emplacing (0): " << bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose() << "\n";
                #endif
            }
            bufferWD.Clear(winSize);
        }

        #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "    newSeed = {" << newSeed.Verbose() << "}\n";
            std::cerr << "    State after: winBuffPos = " << bufferWD.winBuffPos
                            << ", numBasesIn[space] = " << bufferWD.numBasesIn[space]
                            << ", winBuffMinPos = " << bufferWD.winBuffMinPos
                            << ", minimizer = [" << (bufferWD.winBuffMinPos < 0 ? "(no_minimizer)" : bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose()) << "]"
                            << "\n";

        #endif

        // The first time the buffer is filled, find and add the previous
        // distinct minimizer matches.
        if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos) && IsWindowFullFirstTime(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize)) {
            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    First window buffer full!\n";
                std::cerr << "      - minimizer = [" << (bufferWD.winBuffMinPos < 0 ? "(no_minimizer)" : bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose()) << "\n";
                std::cerr << "      - First circular part:\n";
            #endif

            const auto& minSeed = bufferWD.winBuff[bufferWD.winBuffMinPos];
            // First portion of the circular buffer.
            for (int32_t j = (bufferWD.winBuffPos + 1); j < winSize; ++j) {
                const bool isAlt = IsAlternativeMinimizer(minSeed, bufferWD.winBuff[j]);

                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "        - [j = " << j << "] isAlt = " << isAlt << "; bufferWD.winBuff[j] = {" << bufferWD.winBuff[j] << "}, minSeed = {" << minSeed << "}\n";
                #endif

                if (isAlt) {
                    minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                    #ifdef DEBUG_GENERATE_MINIMIZERS
                        std::cerr << "          - Emplacing (1): " << bufferWD.winBuff[j].Verbose() << "\n";
                    #endif
                }
            }

            #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "      - Second circular part:\n";
            #endif

            // Second portion of the circular buffer.
            for (int32_t j = 0; j < bufferWD.winBuffPos; ++j) {
                const bool isAlt = IsAlternativeMinimizer(minSeed, bufferWD.winBuff[j]);

                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "        - [j = " << j << "] isAlt = " << isAlt << "; bufferWD.winBuff[j] = {" << bufferWD.winBuff[j] << "}, minSeed = {" << minSeed << "}\n";
                #endif

                if (isAlt) {
                    minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                    #ifdef DEBUG_GENERATE_MINIMIZERS
                        std::cerr << "          - Emplacing (2): " << bufferWD.winBuff[j].Verbose() << "\n";
                    #endif
                }
            }

            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    Done searching for minimizers in the first window.\n";
            #endif
        }

        /*
         * The first window is always special in this formulation of the minimizer generator:
         * - Any equivalent minimial seed key is written out only after the first window is fully processed.
         *   For any other window, previous minimizers are written out as soon as the new minimizers are found.
         * - This is necessary because until the first window is fully processed, we cannot know that any local minimum
         *   seed key will be the global minimum within that window. This answer will only be known when the first `w`
         *   seeds are written (at the last position of the minimizer window), at which point we need to circle around
         *   the buffer to find all equally minimial keys.
         * - After the first window, we can simply dump the previous minimizer, because we know it was the actual
         *   global minimum of the first window (because the window slides down by 1 base).
        */

        if (newSeed.Valid() && bufferWD.winBuffMinPos < 0) {
            /*
             * No minimum has been set yet, and we found the first valid key.
             * Set the current minimum buffer pos.
             */

            bufferWD.winBuffMinPos = bufferWD.winBuffPos;

            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    (if 1) No minimizers set yet. Setting (1): bufferWD.winBuffMinPos = " << bufferWD.winBuffMinPos << "\n";
            #endif

        } else if (bufferWD.winBuffMinPos >= 0 && newSeed.Valid() && newSeed.key <= bufferWD.winBuff[bufferWD.winBuffMinPos].key) {
            /*
             * We found a new minimum. Write out the previous minimum.
             * In this case, even if we encountered the same minimal key, we will write
             * out the previous occurence of this key, and then set the minimum to the current
             * one. This ensures that all equal minimizers in a window are written.
            */
            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    (else if 2) New minimizer found. Flushing out the old one.\n";
            #endif

            // Special case for handling the first window where the last kmer has an equivalent minimizer
            // seed somewhere before it (e.g. CTCTCT... dinuc repeat). Check out this test: GenerateMinimizers::Short_Dinuc.
            // This is distinct from the `else if` a few lines below because we need to check the actual position and that
            // the seeds have an equivalent key. For the first window, if the key is not equivalent (meaning: it's smaller),
            // then we need to ignore the previous minimum. If it's equivalent, we need to write it to not lose seeds in the
            // first window.
            // For any other window, we write out the seed regardles of if it's equivalent or smaller than before.
            if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos) && IsWindowFullFirstTime(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize) &&
                newSeed.key == bufferWD.winBuff[bufferWD.winBuffMinPos].key) {
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "    - Emplacing (3.1): " << bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose() << "\n";
                #endif
            } else if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos) && IsWindowFull(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize, 0)) {
                /*
                 * Offset in IsWindowFull here is correctly set to 0. This is done in order to check every window
                 * after the first one, since the first window is special. Check the above comment for the explanation.
                */
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "    - Emplacing (3.2): " << bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose() << "\n";
                #endif
            }
            bufferWD.winBuffMinPos = bufferWD.winBuffPos;
            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    Setting (2): bufferWD.winBuffMinPos = " << bufferWD.winBuffMinPos << "\n";
            #endif

        } else if (bufferWD.winBuffMinPos >= 0 && bufferWD.winBuffPos == bufferWD.winBuffMinPos) {
            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    (else if 3) The entire buffer was circled around. Flush current minimizer and find new ones.\n";
            #endif

            /*
             * The entire window has been circled around to the minimum seed key.
             */

            if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos) && IsWindowFull(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize, -1)) {
                minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
                #ifdef DEBUG_GENERATE_MINIMIZERS
                    std::cerr << "    - Emplacing (4): " << bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose() << "\n";
                #endif
            }
            // Find a new minimizer using the two for loops for the two portions of the circular buffer.
            bufferWD.winBuffMinPos = -1;
            // First portion of the circular buffer to find a new minimizer.
            for (int32_t j = (bufferWD.winBuffPos + 1); j < winSize; ++j) {
                if (bufferWD.winBuff[j].Valid() && (bufferWD.winBuffMinPos < 0 || bufferWD.winBuff[j].key <= bufferWD.winBuff[bufferWD.winBuffMinPos].key)) {
                    bufferWD.winBuffMinPos = j;
                }
            }
            // std::cerr << "(1) Looking for new winBuffMinPos, winBuffMinPos = " << winBuffMinPos << "\n";
            // Second portion of the circular buffer to find a new minimizer.
            for (int32_t j = 0; j < bufferWD.winBuffPos; ++j) {
                if (bufferWD.winBuff[j].Valid() && (bufferWD.winBuffMinPos < 0 || bufferWD.winBuff[j].key <= bufferWD.winBuff[bufferWD.winBuffMinPos].key)) {
                    bufferWD.winBuffMinPos = j;
                }
            }
            // std::cerr << "(2) Looking for new winBuffMinPos, winBuffMinPos = " << winBuffMinPos << "\n";
            // Check if the new minimizer is the smallest one.
            if (newSeed.Valid() && IsWindowFull(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize, -1) && (bufferWD.winBuffMinPos < 0 || newSeed.key <= bufferWD.winBuff[bufferWD.winBuffMinPos].key)) {
                bufferWD.winBuffMinPos = bufferWD.winBuffPos;
            }
            #ifdef DEBUG_GENERATE_MINIMIZERS
                std::cerr << "    Setting (3): bufferWD.winBuffMinPos = " << bufferWD.winBuffMinPos << "\n";
            #endif

            // std::cerr << "(3) Looking for new winBuffMinPos, winBuffMinPos = " << winBuffMinPos << "\n";
            // Find and add keys that are identical to the newly found minimizer.
            if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos) && IsWindowFull(bufferWD.numBasesIn[space], winSize, kmerSize, spacing, minSeedSpanSize, -1)) {
                const auto& minSeed = bufferWD.winBuff[bufferWD.winBuffMinPos];
                // First portion of the circular buffer.
                for (int32_t j = (bufferWD.winBuffPos + 1); j < winSize; ++j) {
                    if (IsAlternativeMinimizer(minSeed, bufferWD.winBuff[j])) {
                        minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                        #ifdef DEBUG_GENERATE_MINIMIZERS
                            std::cerr << "    - Emplacing (5): " << bufferWD.winBuff[j].Verbose() << "\n";
                        #endif
                    }
                }
                // Second portion of the circular buffer.
                for (int32_t j = 0; j < bufferWD.winBuffPos; ++j) {
                    if (IsAlternativeMinimizer(minSeed, bufferWD.winBuff[j])) {
                        minimizers.emplace_back(bufferWD.winBuff[j].To128t());
                        #ifdef DEBUG_GENERATE_MINIMIZERS
                            std::cerr << "    - Emplacing (6): " << bufferWD.winBuff[j].Verbose() << "\n";
                        #endif
                    }
                }
            }
        }
        // Store the new seed to the circular buffer.
        bufferWD.winBuff[bufferWD.winBuffPos] = newSeed;
        ++bufferWD.winBuffPos;
        if (bufferWD.winBuffPos == winSize) {
            bufferWD.winBuffPos = 0;
        }

        #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "\n";
        #endif
    }
    if (IsWindowBufferElementValid(bufferWD.winBuff, bufferWD.winBuffMinPos)) {
        minimizers.emplace_back(bufferWD.winBuff[bufferWD.winBuffMinPos].To128t());
        #ifdef DEBUG_GENERATE_MINIMIZERS
            std::cerr << "    - Emplacing (7): " << bufferWD.winBuff[bufferWD.winBuffMinPos].Verbose() << "\n";
        #endif
    }

    #ifdef DEBUG_GENERATE_MINIMIZERS
        std::cerr << "\n";
    #endif

    return 0;
    // clang-format on
}

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        const std::vector<FastaSequenceCached>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC)
{
    retSeeds.clear();
    for (int32_t recordId = 0; recordId < static_cast<int32_t>(targetSeqs.size()); ++recordId) {
        const auto& record = targetSeqs[recordId];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.data());
        const int32_t seqId = record.Id();
        int32_t seqLen = record.size();
        std::vector<PacBio::Pancake::Int128t> newSeeds;
        int rv = GenerateMinimizers(newSeeds, seq, seqLen, 0, seqId, kmerSize, winSize, spacing,
                                    useReverseComplement, useHPC);
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the target sequence, id = " +
                                     std::to_string(recordId));
        retSeeds.insert(retSeeds.end(), newSeeds.begin(), newSeeds.end());
    }
}

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        const std::vector<std::string>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC)
{
    retSeeds.clear();
    for (int32_t seqId = 0; seqId < static_cast<int32_t>(targetSeqs.size()); ++seqId) {
        const auto& record = targetSeqs[seqId];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.data());
        int32_t seqLen = record.size();
        std::vector<PacBio::Pancake::Int128t> newSeeds;
        int rv = GenerateMinimizers(newSeeds, seq, seqLen, 0, seqId, kmerSize, winSize, spacing,
                                    useReverseComplement, useHPC);
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the target sequence, id = " +
                                     std::to_string(seqId));
        retSeeds.insert(retSeeds.end(), newSeeds.begin(), newSeeds.end());
    }
}

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<FastaSequenceCached>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC)
{
    retSeeds.clear();
    retSequenceLengths.clear();
    retSequenceLengths.reserve(targetSeqs.size());
    for (int32_t recordId = 0; recordId < static_cast<int32_t>(targetSeqs.size()); ++recordId) {
        const auto& record = targetSeqs[recordId];
        const uint8_t* seq = reinterpret_cast<const uint8_t*>(record.data());
        const int32_t seqId = record.Id();
        int32_t seqLen = record.size();
        retSequenceLengths.emplace_back(seqLen);
        std::vector<PacBio::Pancake::Int128t> newSeeds;
        int rv = GenerateMinimizers(newSeeds, seq, seqLen, 0, seqId, kmerSize, winSize, spacing,
                                    useReverseComplement, useHPC);
        if (rv)
            throw std::runtime_error("Generating minimizers failed for the target sequence, id = " +
                                     std::to_string(recordId));
        retSeeds.insert(retSeeds.end(), newSeeds.begin(), newSeeds.end());
    }
}

void GenerateMinimizers(std::vector<PacBio::Pancake::Int128t>& retSeeds,
                        std::vector<int32_t>& retSequenceLengths,
                        const std::vector<std::string>& targetSeqs, const int32_t kmerSize,
                        const int32_t winSize, const int32_t spacing,
                        const bool useReverseComplement, const bool useHPC)
{
    std::vector<FastaSequenceCached> targetSeqsCached;
    for (int32_t i = 0; i < static_cast<int32_t>(targetSeqs.size()); ++i) {
        targetSeqsCached.emplace_back(
            FastaSequenceCached(std::to_string(i), targetSeqs[i].c_str(), targetSeqs[i].size(), i));
    }
    GenerateMinimizers(retSeeds, retSequenceLengths, targetSeqsCached, kmerSize, winSize, spacing,
                       useReverseComplement, useHPC);
}

}  // namespace SeedDB
}  // namespace Pancake
}  // namespace PacBio
