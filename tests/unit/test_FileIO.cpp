// Authors: Ivan Sovic

#include <PancakeTestData.h>
#include <gtest/gtest.h>
#include <cstdint>
#include <fstream>
#include <pancake/util/FileIO.hpp>
#include <string>
#include <tuple>
#include <unordered_set>
#include <vector>

namespace FileIOTests {

TEST(Test_FileIO, LoadLinesToVector)
{
    // Input data which will be written and read.
    std::vector<std::string> inData = {"2", "1", "3", "4", "4", "4", "5"};
    std::vector<std::string> expected = inData;

    // The input file for the parsing function.
    std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/in.data";

    // Write the test input file.
    {
        std::ofstream ofs(tmpInFile);
        for (const auto& line : inData) {
            ofs << line << "\n";
        }
    }

    // Run unit under test.
    std::vector<std::string> results = PacBio::Pancake::LoadLinesToVector(tmpInFile);

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, LoadLinesToSet)
{
    // Input data which will be written and read.
    std::vector<std::string> inData = {"2", "1", "3", "4", "4", "4", "5"};
    std::unordered_set<std::string> expected = {"2", "1", "3", "4", "5"};

    // The input file for the parsing function.
    std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/in.data";

    // Write the test input file.
    {
        std::ofstream ofs(tmpInFile);
        for (const auto& line : inData) {
            ofs << line << "\n";
        }
    }

    // Run unit under test.
    std::unordered_set<std::string> results = PacBio::Pancake::LoadLinesToSet(tmpInFile);

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, LoadLinesToOrderedSet)
{
    // Input data which will be written and read.
    std::vector<std::string> inData = {"2", "1", "3", "4", "4", "4", "5"};
    std::set<std::string> expected = {"2", "1", "3", "4", "5"};

    // The input file for the parsing function.
    std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/in.data";

    // Write the test input file.
    {
        std::ofstream ofs(tmpInFile);
        for (const auto& line : inData) {
            ofs << line << "\n";
        }
    }

    // Run unit under test.
    std::set<std::string> results = PacBio::Pancake::LoadLinesToOrderedSet(tmpInFile);

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, FormatIs)
{
    // Input data.
    // clang-format off
    // Elements of the tuple:
    //  <path, isFasta, isFastq, isBam, isXml, isFofn, isSeqDB>
    std::vector<std::tuple<std::string, bool, bool, bool, bool, bool, bool>> inData = {
        {"test.fasta", true, false, false, false, false, false},
        {"test.fa", true, false, false, false, false, false},
        {"test.fasta.gz", true, false, false, false, false, false},
        {"test.fa.gz", true, false, false, false, false, false},
        {"test.FASTA", true, false, false, false, false, false},
        {"test.fastq", false, true, false, false, false, false},
        {"test.fq", false, true, false, false, false, false},
        {"test.fastq.gz", false, true, false, false, false, false},
        {"test.fq.gz", false, true, false, false, false, false},
        {"test.FASTQ", false, true, false, false, false, false},
        {"test.bam", false, false, true, false, false, false},
        {"test.BAM", false, false, true, false, false, false},
        {"test.xml", false, false, false, true, false, false},
        {"test.XML", false, false, false, true, false, false},
        {"test.fofn", false, false, false, false, true, false},
        {"test.FOFN", false, false, false, false, true, false},
        {"test.seqdb", false, false, false, false, false, true},
        {"test.SeqDB", false, false, false, false, false, true},
        {"something_else", false, false, false, false, false, false},
        {"", false, false, false, false, false, false},
        {"test.", false, false, false, false, false, false},
    };
    // clang-format on

    // Batch test.
    for (const auto& inPair : inData) {
        const auto& inFile = std::get<0>(inPair);
        const auto& isFasta = std::get<1>(inPair);
        const auto& isFastq = std::get<2>(inPair);
        const auto& isBam = std::get<3>(inPair);
        const auto& isXml = std::get<4>(inPair);
        const auto& isFofn = std::get<5>(inPair);
        const auto& isSeqDB = std::get<6>(inPair);

        {
            bool expected = isFasta;
            bool result = PacBio::Pancake::FormatIsFasta(inFile);
            EXPECT_EQ(expected, result);
        }
        {
            bool expected = isFastq;
            bool result = PacBio::Pancake::FormatIsFastq(inFile);
            EXPECT_EQ(expected, result);
        }
        {
            bool expected = isBam;
            bool result = PacBio::Pancake::FormatIsBam(inFile);
            EXPECT_EQ(expected, result);
        }
        {
            bool expected = isXml;
            bool result = PacBio::Pancake::FormatIsXml(inFile);
            EXPECT_EQ(expected, result);
        }
        {
            bool expected = isFofn;
            bool result = PacBio::Pancake::FormatIsFofn(inFile);
            EXPECT_EQ(expected, result);
        }
        {
            bool expected = isSeqDB;
            bool result = PacBio::Pancake::FormatIsSeqDB(inFile);
            EXPECT_EQ(expected, result);
        }
    }
}

TEST(Test_FileIO, ExpandInputFileList_EmptyInput)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {};
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {};

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, true);

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_OnlyNonComposite)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        "test1.bam",
    };
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Bam, "test1.bam"},
    };

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, true);

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_SingleFofnAndOneNoncomposite)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn",
        "test2.bam",
    };

    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Fasta, "test1a.fasta"},
        {PacBio::Pancake::SequenceFormat::Fastq, "test1b.fastq"},
        {PacBio::Pancake::SequenceFormat::Bam, "test2.bam"},
    };
    std::sort(expected.begin(), expected.end());

    // Create the FOFNs.
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << "test1a.fasta"
            << "\n";
        ofs << "test1b.fastq";
    }

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, true);
    std::sort(results.begin(), results.end());

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_CompositeMix)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn",
        PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test2.fofn",
        "test4.bam",
    };

    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Fasta, "test1a.fasta"},
        {PacBio::Pancake::SequenceFormat::Fastq, "test1b.fastq"},
        {PacBio::Pancake::SequenceFormat::SeqDB, "test3a.seqdb"},
        {PacBio::Pancake::SequenceFormat::Fasta, "test3b.fasta"},
        {PacBio::Pancake::SequenceFormat::Fastq, "test3c.fastq"},
        {PacBio::Pancake::SequenceFormat::Bam, "test4.bam"},
    };
    std::sort(expected.begin(), expected.end());

    // Create the FOFNs.
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << "test1a.fasta"
            << "\n";
        ofs << "test1b.fastq";
    }
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test2.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test3.fofn"
            << "\n";
    }
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test3.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << "test3a.seqdb"
            << "\n";
        ofs << "test3b.fasta"
            << "\n";
        ofs << "test3c.fastq"
            << "\n";
    }

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, true);
    std::sort(results.begin(), results.end());

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_CircularFofnsShouldThrow)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn",
    };

    // Create the FOFNs.
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test2.fofn"
            << "\n";
    }
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test2.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn"
            << "\n";
    }

    // Run unit under  test.
    EXPECT_THROW({ PacBio::Pancake::ExpandInputFileList(inFiles, true); }, std::runtime_error);
}

TEST(Test_FileIO, ExpandInputFileList_FofnWithAnXmlAndExpansionToBAM)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn",
    };

    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Bam,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreads1.bam"},
        {PacBio::Pancake::SequenceFormat::Bam,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreads2.bam"},
        {PacBio::Pancake::SequenceFormat::Bam,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreads3.bam"},
    };
    std::sort(expected.begin(), expected.end());

    // Create the FOFNs.
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreadset.xml";
    }

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, true);
    std::sort(results.begin(), results.end());

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_NoFofnJustXmlWithExpansionToBAM)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreadset.xml",
    };

    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Bam,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreads1.bam"},
        {PacBio::Pancake::SequenceFormat::Bam,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreads2.bam"},
        {PacBio::Pancake::SequenceFormat::Bam,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreads3.bam"},
    };
    std::sort(expected.begin(), expected.end());

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, true);
    std::sort(results.begin(), results.end());

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_FofnWithAnXmlAndNoExpansionToBAM)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn",
    };

    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Xml,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreadset.xml"},
    };
    std::sort(expected.begin(), expected.end());

    // Create the FOFNs.
    {
        // The input file for the parsing function.
        std::string tmpInFile = PacBio::PancakeTestsConfig::GeneratedData_Dir + "/test1.fofn";
        std::ofstream ofs(tmpInFile);
        ofs << PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreadset.xml";
    }

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, false);
    std::sort(results.begin(), results.end());

    // Verify.
    EXPECT_EQ(expected, results);
}

TEST(Test_FileIO, ExpandInputFileList_NoFofnJustXmlWithNoExpansionToBAM)
{
    // Input data which will be written and read.
    std::vector<std::string> inFiles = {
        PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreadset.xml",
    };

    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> expected = {
        {PacBio::Pancake::SequenceFormat::Xml,
         PacBio::PancakeTestsConfig::Data_Dir + "/seqdb-writer/bam/subreadset.xml"},
    };
    std::sort(expected.begin(), expected.end());

    // Run unit under  test.
    std::vector<std::pair<PacBio::Pancake::SequenceFormat, std::string>> results =
        PacBio::Pancake::ExpandInputFileList(inFiles, false);
    std::sort(results.begin(), results.end());

    // Verify.
    EXPECT_EQ(expected, results);
}
}  // namespace FileIOTests
