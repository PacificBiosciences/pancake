// Authors: Ivan Sovic

#include "SeqFetchWorkflow.h"
#include <pacbio/util/FileIO.h>
#include "SeqFetchSettings.h"

#include <algorithm>
#include <cinttypes>
#include <iostream>
#include <set>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

#include <pacbio/pancake/SeqDBIndexCache.h>
#include <pacbio/pancake/SeqDBReader.h>
#include <pacbio/util/RunLengthEncoding.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/IndexedFastqReader.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>
#include <boost/algorithm/string/predicate.hpp>

#include <pbcopper/logging/LogLevel.h>
#include <pbcopper/logging/Logging.h>

namespace PacBio {
namespace Pancake {

void WriteSeq(std::ostream& os, const std::string& seqName, const std::string& seq,
              const std::string& quals, char dummyQV, bool useProvidedQuals,
              const SeqFetchOutFormat& outFmt)
{
    if (outFmt == SeqFetchOutFormat::Fastq) {
        os << "@" << seqName << "\n" << seq << "\n+\n";
        if (useProvidedQuals) {
            os << quals << "\n";
        } else {
            os << std::string(seq.size(), dummyQV) << "\n";
        }
    } else if (outFmt == SeqFetchOutFormat::Fasta) {
        os << ">" << seqName << "\n" << seq << "\n";
    } else {
        throw std::runtime_error("Unknown output format!");
    }
}

void WriteRLE(std::ostream& os, const std::string& seqName, std::vector<int32_t> hpcToSeqCoords)
{
    os << ">" << seqName << "\n";
    if (hpcToSeqCoords.size() > 0) {
        os << hpcToSeqCoords[0];
    }
    for (size_t i = 1; i < hpcToSeqCoords.size(); ++i) {
        os << ",";
        os << hpcToSeqCoords[i];
    }
    os << "\n";
}

void WriteSeqAndRLE(std::shared_ptr<std::ostream>& osPtr, std::shared_ptr<std::ostream>& osRlePtr,
                    const std::string& seqName, const std::string& seq, const std::string& quals,
                    char dummyQV, bool useProvidedQuals, const SeqFetchOutFormat& outFmt,
                    bool useHPC, bool useRLE)
{
    std::string seqRLE;
    std::vector<int32_t> seqToHPCCoords;
    std::vector<int32_t> hpcToSeqCoords;

    // Compute the run length encoding when either useHPC or useRLE are specified.
    // In case of useRLE, only the coordinates will be used, and in case of
    // the useHPC only the sequence will be used.
    if (useHPC || useRLE) {
        RunLengthEncoding(seq, seqRLE, seqToHPCCoords, hpcToSeqCoords);
    }

    // Only use the HPC sequence if useHPC is specified.
    const std::string& seqFinal = useHPC ? seqRLE : seq;

    // Write the sequence and quals.
    WriteSeq(*osPtr, seqName, seqFinal, quals, dummyQV, useProvidedQuals, outFmt);

    // Write the RLE if required.
    if (useRLE) {
        WriteRLE(*osRlePtr, seqName, hpcToSeqCoords);
    }
}

void FetchFromFasta(std::shared_ptr<std::ostream>& osPtr, std::shared_ptr<std::ostream>& osRlePtr,
                    std::vector<std::string>& foundSeqs, const std::string& inFile,
                    const std::set<std::string>& remainingToFind, const char dummyQV,
                    const PacBio::Pancake::SeqFetchOutFormat& outFormat, bool useHPC, bool useRLE)
{
    BAM::IndexedFastaReader reader{inFile};
    PBLOG_INFO << "Num sequences in file: " << reader.NumSequences();
    if (reader.NumSequences() == 0) {
        throw std::runtime_error("Index does not exist for file '" + inFile + "'!");
    }
    for (const auto& seqName : remainingToFind) {
        if (!reader.HasSequence(seqName)) continue;
        std::string seq = reader.Subsequence(seqName.c_str());
        std::string qual;
        WriteSeqAndRLE(osPtr, osRlePtr, seqName, seq, qual, dummyQV, false, outFormat, useHPC,
                       useRLE);
        foundSeqs.emplace_back(seqName);
    }
}

void FetchFromFastq(std::shared_ptr<std::ostream>& osPtr, std::shared_ptr<std::ostream>& osRlePtr,
                    std::vector<std::string>& foundSeqs, const std::string& inFile,
                    const std::set<std::string>& remainingToFind, const char dummyQV,
                    const PacBio::Pancake::SeqFetchOutFormat& outFormat, bool useHPC, bool useRLE)
{
    BAM::IndexedFastqReader reader{inFile};
    PBLOG_INFO << "Num sequences in file: " << reader.NumSequences();
    if (reader.NumSequences() == 0) {
        throw std::runtime_error("Index does not exist for file '" + inFile + "'!");
    }
    for (const auto& seqName : remainingToFind) {
        if (!reader.HasSequence(seqName)) continue;
        int32_t seqLen = reader.SequenceLength(seqName);
        std::pair<std::string, PacBio::BAM::QualityValues> seqQualPair =
            reader.Subsequence(seqName, 0, seqLen);
        const auto& seq = std::get<0>(seqQualPair);
        std::string qual = std::get<1>(seqQualPair).Fastq();
        foundSeqs.emplace_back(seqName);
        WriteSeqAndRLE(osPtr, osRlePtr, seqName, seq, qual, dummyQV, true, outFormat, useHPC,
                       useRLE);
    }
}

void FetchFromBam(std::shared_ptr<std::ostream>& osPtr, std::shared_ptr<std::ostream>& osRlePtr,
                  std::vector<std::string>& foundSeqs, const std::string& inFile,
                  const std::set<std::string>& remainingToFind, const char dummyQV,
                  const PacBio::Pancake::SeqFetchOutFormat& outFormat, bool useHPC, bool useRLE)
{
    PacBio::BAM::BamFile bamFile(inFile);
    bamFile.EnsurePacBioIndexExists();

    std::vector<std::string> names;
    names.insert(names.end(), remainingToFind.begin(), remainingToFind.end());
    PacBio::BAM::PbiQueryNameFilter filter{names};
    PacBio::BAM::PbiFilterQuery query{filter, bamFile};
    for (auto record : query) {
        auto quals = record.Qualities().Fastq();
        WriteSeqAndRLE(osPtr, osRlePtr, record.FullName(), record.Sequence(), quals, dummyQV,
                       quals.size() > 0, outFormat, useHPC, useRLE);
        foundSeqs.emplace_back(record.FullName());
    }
}

void FetchFromSeqDB(std::shared_ptr<std::ostream>& osPtr, std::shared_ptr<std::ostream>& osRlePtr,
                    std::vector<std::string>& foundSeqs, const std::string& inFile,
                    const std::set<std::string>& remainingToFind, const char dummyQV,
                    const bool writeIds, const PacBio::Pancake::SeqFetchOutFormat& outFormat,
                    bool useHPC, bool useRLE)
{
    std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
        PacBio::Pancake::LoadSeqDBIndexCache(inFile);

    PacBio::Pancake::SeqDBReader reader(seqDBCache);

    PacBio::Pancake::FastaSequenceId record;
    for (const auto& seqName : remainingToFind) {
        bool rv = reader.GetSequence(record, seqName);
        if (rv == false) {
            continue;
        }
        foundSeqs.emplace_back(seqName);
        if (writeIds) {
            char buff[50];
            sprintf(buff, "%09" PRIi64, record.Id());
            WriteSeqAndRLE(osPtr, osRlePtr, std::string(buff), record.Bases(), std::string(),
                           dummyQV, false, outFormat, useHPC, useRLE);
        } else {
            WriteSeqAndRLE(osPtr, osRlePtr, seqName, record.Bases(), std::string(), dummyQV, false,
                           outFormat, useHPC, useRLE);
        }
    }
}

int SeqFetchWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqFetchSettings settings{options};

    // Expand FOFNs and determine the formats of input files.
    std::vector<std::pair<SequenceFormat, std::string>> inFiles =
        ExpandInputFileList(settings.InputFiles, false);

    // Parse the specified list of sequences to be fetched.
    std::set<std::string> seqNamesToFind = LoadLinesToOrderedSet(settings.InputFetchListFile);

    // Verify that the settings.WriteIds is used properly.
    for (const auto& filePair : inFiles) {
        const auto& inFmt = std::get<0>(filePair);
        const auto& inFile = std::get<1>(filePair);

        if (settings.WriteIds && inFmt != SequenceFormat::SeqDB) {
            std::ostringstream oss;
            oss << "Cannot use the --write-ids option with input files which are not in the SeqDB "
                   "format. Offending file: '"
                << inFile << "'.";
            throw std::runtime_error(oss.str());
        }
    }

    // If the input list is composed of sequence IDs, then the alias SeqDB
    // can be used to get the actual sequence names corresponding to those IDs.
    if (!settings.AliasSeqDBFile.empty()) {
        std::shared_ptr<PacBio::Pancake::SeqDBIndexCache> seqDBCache =
            PacBio::Pancake::LoadSeqDBIndexCache(settings.AliasSeqDBFile);

        // The SeqDB is defined such that the ID matches the ordinal number of the seq line.
        std::vector<std::string> actualNames;
        for (const auto& idStr : seqNamesToFind) {
            int32_t id = std::stoi(idStr);
            const auto& sl = seqDBCache->GetSeqLine(id);
            actualNames.emplace_back(sl.header);
        }

        // Update the name list.
        seqNamesToFind.clear();
        seqNamesToFind.insert(actualNames.begin(), actualNames.end());
    }

    // Output file stream. If the file name is "-", then write to stdout.
    std::shared_ptr<std::ostream> osPtr(&std::cout, [](void*) {});
    if (settings.OutputFile.size() > 0 && settings.OutputFile != "-") {
        osPtr = std::shared_ptr<std::ostream>(new std::ofstream(settings.OutputFile));
        PBLOG_INFO << "Output is to file: " << settings.OutputFile;
    } else {
        if (settings.UseRLE) {
            throw std::runtime_error(
                "Cannot output to sequences to stdout and write a .rle file. Please specify a "
                "concrete output file.");
        }
        PBLOG_INFO << "Output is to stdout.";
    }

    std::shared_ptr<std::ostream> osRlePtr = nullptr;
    if (settings.UseRLE) {
        osRlePtr = std::shared_ptr<std::ostream>(new std::ofstream(settings.OutputFile + ".rle"));
    }

    PBLOG_INFO << "Starting to fetch.";
    auto remainingToFind = seqNamesToFind;
    for (const auto& filePair : inFiles) {
        if (remainingToFind.empty()) {
            break;
        }

        const auto& inFmt = std::get<0>(filePair);
        const auto& inFile = std::get<1>(filePair);

        PBLOG_INFO << "Looking in file: " << inFile;

        std::vector<std::string> foundSeqs;

        if (inFmt == SequenceFormat::Fasta) {
            FetchFromFasta(osPtr, osRlePtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                           settings.OutputFormat, settings.UseHPC, settings.UseRLE);

        } else if (inFmt == SequenceFormat::Fastq) {
            FetchFromFastq(osPtr, osRlePtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                           settings.OutputFormat, settings.UseHPC, settings.UseRLE);

        } else if (inFmt == SequenceFormat::SeqDB) {
            FetchFromSeqDB(osPtr, osRlePtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                           settings.WriteIds, settings.OutputFormat, settings.UseHPC,
                           settings.UseRLE);

        } else if (inFmt == SequenceFormat::Bam) {
            FetchFromBam(osPtr, osRlePtr, foundSeqs, inFile, remainingToFind, settings.DummyQV,
                         settings.OutputFormat, settings.UseHPC, settings.UseRLE);
        }

        // Remove the found sequences from the remaining set.
        for (const auto& found : foundSeqs) {
            remainingToFind.erase(found);
        }
    }

    PBLOG_INFO << "Done!";
    PBLOG_INFO << "Found sequences: "
               << (static_cast<int32_t>(seqNamesToFind.size()) -
                   static_cast<int32_t>(remainingToFind.size()))
               << " / " << seqNamesToFind.size() << " .";

    if (settings.FailOnMissingQueries && remainingToFind.size() > 0) {
        std::ostringstream oss;
        oss << "Not all queries were found in the provided input files! Found sequences: "
            << (static_cast<int32_t>(seqNamesToFind.size()) -
                static_cast<int32_t>(remainingToFind.size()))
            << " / " << seqNamesToFind.size() << " .";
        throw std::runtime_error(oss.str());
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
