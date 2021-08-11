// Author: Ivan Sovic

#ifndef PANCAKE_FILE_IO_H
#define PANCAKE_FILE_IO_H

#include <cstdint>
#include <set>
#include <string>
#include <unordered_set>
#include <vector>

namespace PacBio {
namespace Pancake {

/// \brief Enum class with all supported input sequence formats.
enum class SequenceFormat
{
    Fasta,
    Fastq,
    Bam,
    Xml,
    Fofn,
    SeqDB,
    Unknown,
};

/// \brief Parses a file line by line, and returns a vector of strings, one per line.
std::vector<std::string> LoadLinesToVector(const std::string& listPath);

/// \brief Parses a file line by line, and stores each line into a set.
///         Useful for parsing sets of reads or similar identifiers.
std::unordered_set<std::string> LoadLinesToSet(const std::string& listPath);

/// \brief Parses a file line by line, and stores each line into a set.
///         Useful for parsing sets of reads or similar identifiers.
std::set<std::string> LoadLinesToOrderedSet(const std::string& listPath);

/// \brief Checks if the provided file path has a FASTA/gzipped FASTA type of extension.
bool FormatIsFasta(const std::string& fn);

/// \brief Checks if the provided file path has a FASTQ/gzipped FASTQ type of extension.
bool FormatIsFastq(const std::string& fn);

/// \brief Checks if the provided file path has a FOFN extension.
bool FormatIsFofn(const std::string& fn);

/// \brief Checks if the provided file path has a BAM extension.
bool FormatIsBam(const std::string& fn);

/// \brief Checks if the provided file path has an XML extension.
bool FormatIsXml(const std::string& fn);

/// \brief Checks if the provided file path has a SeqDBe xtension.
bool FormatIsSeqDB(const std::string& fn);

/// \brief Parses the format of a file based on it's extension.
SequenceFormat ParseFormat(const std::string& filename);

/// \brief Converts the sequence format to a std::string representation.
std::string SequenceFormatToString(const SequenceFormat fmt);

/// \brief Recursively expands all composite files in the list (FOFN, XML) until
///         there are no more composite files left.
/// \returns Vector of pairs of final non-composite files and their file formats.
std::vector<std::pair<SequenceFormat, std::string>> ExpandInputFileList(
    const std::vector<std::string>& inFiles, bool expandXML);

}  // namespace Pancake
}  // namespace PacBio

#endif  // PANCAKE_FILE_IO_H
