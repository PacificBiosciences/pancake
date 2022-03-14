// Authors: Ivan Sovic

#include "SeqDBWorkflow.h"
#include "SeqDBSettings.h"

#include <pancake/SeqDBWriter.hpp>
#include <pancake/util/FileIO.hpp>

#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/FastaReader.h>
#include <pbbam/FastqReader.h>
#include <pbbam/PbiIndexedBamReader.h>
#include <boost/algorithm/string/predicate.hpp>

#include <iostream>
#include <string>
#include <vector>

namespace PacBio {
namespace Pancake {

int SeqDBWorkflow::Runner(const PacBio::CLI_v2::Results& options)
{
    SeqDBSettings settings{options};

    auto writer = PacBio::Pancake::CreateSeqDBWriter(settings.OutputPrefix,
                                                     settings.CompressionLevel, settings.BufferSize,
                                                     settings.BlockSize, settings.SplitBlocks);

    // Expand FOFNs and determine the formats of input files.
    std::vector<std::pair<SequenceFormat, std::string>> inputFiles =
        ExpandInputFileList(settings.InputFiles, false);

    for (const auto& inFilePair : inputFiles) {
        const auto& inFmt = inFilePair.first;
        const auto& inFile = inFilePair.second;
        if (inFmt == SequenceFormat::Fasta) {
            BAM::FastaReader inReader{inFile};
            BAM::FastaSequence record;
            while (inReader.GetNext(record)) {
                writer->AddSequence(record.Name(), record.Bases());
            }
        } else if (inFmt == SequenceFormat::Fastq) {
            BAM::FastqReader inReader{inFile};
            BAM::FastqSequence record;
            while (inReader.GetNext(record)) {
                writer->AddSequence(record.Name(), record.Bases());
            }
        } else if (inFmt == SequenceFormat::Bam) {
            BAM::BamReader inputBamReader{inFile};
            for (const auto& bam : inputBamReader)
                writer->AddSequence(bam.FullName(), bam.Sequence());
        } else if (inFmt == SequenceFormat::Xml) {
            BAM::DataSet dataset{inFile};
            const PacBio::BAM::PbiIndexCache pbiCache = PacBio::BAM::MakePbiIndexCache(dataset);
            const PacBio::BAM::PbiFilter filter = PacBio::BAM::PbiFilter::FromDataSet(dataset);
            size_t fileId = 0;
            for (const PacBio::BAM::BamFile& bam : dataset.BamFiles()) {
                const std::shared_ptr<PacBio::BAM::PbiRawData>& pbiIndex = pbiCache->at(fileId);
                PacBio::BAM::PbiIndexedBamReader reader{filter, bam, pbiIndex};
                for (const auto& record : reader) {
                    writer->AddSequence(record.FullName(), record.Sequence());
                }
                ++fileId;
            }
        } else {
            throw std::runtime_error("Unknown input file extension for file: '" + inFile + "'.");
        }
    }

    return EXIT_SUCCESS;
}

}  // namespace Pancake
}  // namespace PacBio
