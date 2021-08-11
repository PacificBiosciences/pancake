// Authors: Ivan Sovic

#include <pacbio/util/FileIO.h>
#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/IndexedFastaReader.h>
#include <pbbam/IndexedFastqReader.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiFilterTypes.h>
#include <boost/algorithm/string/predicate.hpp>
#include <deque>
#include <fstream>

namespace PacBio {
namespace Pancake {

bool FormatIsFasta(const std::string& fn)
{
    return boost::algorithm::iends_with(fn, ".fasta") ||
           boost::algorithm::iends_with(fn, ".fasta.gz") ||
           boost::algorithm::iends_with(fn, ".fa") || boost::algorithm::iends_with(fn, ".fa.gz");
}

bool FormatIsFastq(const std::string& fn)
{
    return boost::algorithm::iends_with(fn, ".fastq") ||
           boost::algorithm::iends_with(fn, ".fastq.gz") ||
           boost::algorithm::iends_with(fn, ".fq") || boost::algorithm::iends_with(fn, ".fq.gz");
}

bool FormatIsFofn(const std::string& fn) { return boost::algorithm::iends_with(fn, ".fofn"); }

bool FormatIsBam(const std::string& fn) { return boost::algorithm::iends_with(fn, ".bam"); }

bool FormatIsXml(const std::string& fn) { return boost::algorithm::iends_with(fn, ".xml"); }

bool FormatIsSeqDB(const std::string& fn) { return boost::algorithm::iends_with(fn, ".seqdb"); }

SequenceFormat ParseFormat(const std::string& filename)
{
    if (FormatIsFasta(filename)) {
        return SequenceFormat::Fasta;
    } else if (FormatIsFastq(filename)) {
        return SequenceFormat::Fastq;
    } else if (FormatIsBam(filename)) {
        return SequenceFormat::Bam;
    } else if (FormatIsXml(filename)) {
        return SequenceFormat::Xml;
    } else if (FormatIsFofn(filename)) {
        return SequenceFormat::Fofn;
    } else if (FormatIsSeqDB(filename)) {
        return SequenceFormat::SeqDB;
    }
    return SequenceFormat::Unknown;
}

std::string SequenceFormatToString(const SequenceFormat fmt)
{
    if (fmt == SequenceFormat::Fasta) {
        return "fasta";
    } else if (fmt == SequenceFormat::Fastq) {
        return "fastq";
    } else if (fmt == SequenceFormat::Bam) {
        return "bam";
    } else if (fmt == SequenceFormat::Xml) {
        return "xml";
    } else if (fmt == SequenceFormat::Fofn) {
        return "fofn";
    } else if (fmt == SequenceFormat::SeqDB) {
        return "seqdb";
    }
    return "unknown";
}

std::vector<std::string> LoadLinesToVector(const std::string& listPath)
{
    std::vector<std::string> ret;
    std::ifstream ifs(listPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open file '" + listPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        ret.emplace_back(line);
    }
    return ret;
}

std::unordered_set<std::string> LoadLinesToSet(const std::string& listPath)
{
    std::unordered_set<std::string> ret;
    std::ifstream ifs(listPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open file '" + listPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        ret.emplace(line);
    }
    return ret;
}

std::set<std::string> LoadLinesToOrderedSet(const std::string& listPath)
{
    std::set<std::string> ret;
    std::ifstream ifs(listPath);
    if (ifs.is_open() == false) {
        throw std::runtime_error("Could not open file '" + listPath + "'!");
    }
    std::string line;
    while (std::getline(ifs, line)) {
        ret.emplace(line);
    }
    return ret;
}

std::vector<std::pair<SequenceFormat, std::string>> ExpandInputFileList(
    const std::vector<std::string>& inFiles, bool expandXML)
{
    std::deque<std::string> que;

    for (const auto& file : inFiles) {
        que.emplace_back(file);
    }

    std::vector<std::pair<SequenceFormat, std::string>> retFiles;
    std::unordered_set<std::string> visited;

    while (que.size() > 0) {
        auto inFile = que.front();
        que.pop_front();

        auto fmt = ParseFormat(inFile);

        if (fmt == SequenceFormat::Fofn) {
            if (visited.find(inFile) != visited.end()) {
                throw std::runtime_error("Circular FOFN definition detected! File '" + inFile +
                                         "' was already parsed.");
            }

            std::vector<std::string> files = LoadLinesToVector(inFile);
            for (const auto& newFile : files) {
                que.emplace_back(newFile);
            }

        } else if (fmt == SequenceFormat::Xml && expandXML) {
            if (visited.find(inFile) != visited.end()) {
                throw std::runtime_error("Circular XML definition detected! File '" + inFile +
                                         "' was already parsed.");
            }

            BAM::DataSet dataset{inFile};
            const auto& bamFiles = dataset.BamFiles();
            for (const auto& bamFile : bamFiles) {
                auto newFile = bamFile.Filename();
                que.emplace_back(newFile);
            }

        } else {
            retFiles.emplace_back(std::make_pair(fmt, inFile));
        }

        visited.emplace(inFile);
    }

    return retFiles;
}

}  // namespace Pancake
}  // namespace PacBio
