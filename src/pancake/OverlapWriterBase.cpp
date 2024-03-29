// Authors: Ivan Sovic

#include <pancake/OverlapWriterBase.hpp>

namespace PacBio {
namespace Pancake {

void OverlapWriterBase::PrintOverlapAsIPAOvl(std::FILE* fpOut, const Overlap& ovl,
                                             const std::string& Aname, const std::string& Bname,
                                             bool writeIds, bool writeCigar)
{
    /**
     * Aid Bid score idt Arev Astart Aend Alen Brev Bstart Bend Blen Atype Btype in_phase cigar Avars Bvars label
     */

    double identity = static_cast<double>(ovl.Identity);
    if (identity == 0.0 && ovl.EditDistance >= 0.0) {
        const double editDist = ovl.EditDistance;
        const double qSpan = ovl.ASpan();
        const double tSpan = ovl.BSpan();
        const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
        const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
        identity = std::min(identityQ, identityT);
    }

    // The format specifies coordinates always in the FWD strand.
    const int32_t tStart = ovl.BstartFwd();
    const int32_t tEnd = ovl.BendFwd();
    const int32_t tIsRev = ovl.Brev;
    const int32_t tLen = ovl.Blen;
    const std::string AtypeStr = OverlapTypeToStringSingleChar(ovl.Atype);
    const std::string BtypeStr = OverlapTypeToStringSingleChar(ovl.Btype);

    // [1-14] First 12 columns are the same as in M4 + add the overlap type info for A and B reads.
    if (writeIds) {
        std::fprintf(fpOut, "%09d %09d", ovl.Aid, ovl.Bid);
    } else {
        std::fprintf(fpOut, "%s %s", Aname.c_str(), Bname.c_str());
    }
    std::fprintf(fpOut, " %d %.4lf %d %d %d %d %d %d %d %d %s %s", static_cast<int32_t>(ovl.Score),
                 100.0 * identity, static_cast<int32_t>(ovl.Arev), ovl.Astart, ovl.Aend, ovl.Alen,
                 static_cast<int32_t>(tIsRev), tStart, tEnd, tLen, AtypeStr.c_str(),
                 BtypeStr.c_str());

    // [15] In-phase column.
    std::fprintf(fpOut, " u");

    // [16] Write the CIGAR only if specified, for speed.
    if (writeCigar) {
        if (ovl.Cigar.empty()) {
            std::fprintf(fpOut, " *");
        } else {
            std::fprintf(fpOut, " ");
            for (const auto& op : ovl.Cigar) {
                std::fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
            }
        }

        // [17] Write the A-read variant string, in the fwd orientation of the A-read.
        // Variant string is a list of variant bases for every non-match CIGAR operation.
        if (ovl.Avars.empty()) {
            std::fprintf(fpOut, " *");
        } else {
            if (ovl.Arev) {
                const auto vars = Pancake::ReverseComplement(ovl.Avars, 0, ovl.Avars.size());
                std::fprintf(fpOut, " %s", vars.c_str());
            } else {
                std::fprintf(fpOut, " %s", ovl.Avars.c_str());
            }
        }

        // [18] Write the B-read variant string, in the fwd orientation of the B-read.
        if (ovl.Bvars.empty()) {
            std::fprintf(fpOut, " *");
        } else {
            if (ovl.Brev) {
                const auto vars = Pancake::ReverseComplement(ovl.Bvars, 0, ovl.Bvars.size());
                std::fprintf(fpOut, " %s", vars.c_str());
            } else {
                std::fprintf(fpOut, " %s", ovl.Bvars.c_str());
            }
        }

    } else {
        std::fprintf(fpOut, " * * *");
    }

    // [19] Write the additional label column which is unused for now.
    std::fprintf(fpOut, " *");

    std::fprintf(fpOut, "\n");
}

void OverlapWriterBase::PrintOverlapAsM4(std::FILE* fpOut, const Overlap& ovl,
                                         const std::string& Aname, const std::string& Bname,
                                         bool writeIds, bool writeCigar)
{
    double identity = static_cast<double>(ovl.Identity);
    if (identity == 0.0 && ovl.EditDistance >= 0.0) {
        const double editDist = ovl.EditDistance;
        const double qSpan = ovl.ASpan();
        const double tSpan = ovl.BSpan();
        const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
        const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
        identity = std::min(identityQ, identityT);
    }

    // The format specifies coordinates always in the FWD strand.
    const int32_t tStart = ovl.BstartFwd();
    const int32_t tEnd = ovl.BendFwd();
    const int32_t tIsRev = ovl.Brev;
    const int32_t tLen = ovl.Blen;
    const std::string AtypeStr = OverlapTypeToString(ovl.Atype);

    if (writeIds) {
        std::fprintf(fpOut, "%09d %09d", ovl.Aid, ovl.Bid);
    } else {
        std::fprintf(fpOut, "%s %s", Aname.c_str(), Bname.c_str());
    }

    std::fprintf(fpOut, " %d %.2lf %d %d %d %d %d %d %d %d %s", static_cast<int32_t>(ovl.Score),
                 100.0 * identity, static_cast<int32_t>(ovl.Arev), ovl.Astart, ovl.Aend, ovl.Alen,
                 static_cast<int32_t>(tIsRev), tStart, tEnd, tLen, AtypeStr.c_str());

    if (writeCigar) {
        if (ovl.Cigar.empty()) {
            std::fprintf(fpOut, " *");
        } else {
            std::fprintf(fpOut, " ");
            for (const auto& op : ovl.Cigar) {
                std::fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
            }
        }
    }

    std::fprintf(fpOut, "\n");
}

void OverlapWriterBase::PrintOverlapAsPAF(std::FILE* fpOut, const Overlap& ovl,
                                          const std::string& Aname, const std::string& Bname,
                                          bool writeIds, bool writeCigar)
{
    double identity = static_cast<double>(ovl.Identity);
    if (identity == 0.0 && ovl.EditDistance >= 0.0) {
        const double editDist = ovl.EditDistance;
        const double qSpan = ovl.ASpan();
        const double tSpan = ovl.BSpan();
        const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
        const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
        identity = std::min(identityQ, identityT);
    }

    // The format specifies coordinates always in the FWD strand.
    const int32_t tStart = ovl.BstartFwd();
    const int32_t tEnd = ovl.BendFwd();
    const int32_t tIsRev = ovl.Brev;
    const std::string AtypeStr = OverlapTypeToStringSingleChar(ovl.Atype);
    const std::string BtypeStr = OverlapTypeToStringSingleChar(ovl.Btype);
    const int32_t mapq = 60;

    if (writeIds) {
        std::fprintf(fpOut, "%09d\t%d\t%d\t%d\t%c\t%09d\t%d\t%d\t%d\t%d\t%d\t%d", ovl.Aid, ovl.Alen,
                     ovl.Astart, ovl.Aend, (tIsRev ? '-' : '+'), ovl.Bid, ovl.Blen, tStart, tEnd,
                     ovl.ASpan(), ovl.BSpan(), mapq);
    } else {
        std::fprintf(fpOut, "%s\t%d\t%d\t%d\t%c\t%s\t%d\t%d\t%d\t%d\t%d\t%d", Aname.c_str(),
                     ovl.Alen, ovl.Astart, ovl.Aend, (tIsRev ? '-' : '+'), Bname.c_str(), ovl.Blen,
                     tStart, tEnd, ovl.ASpan(), ovl.BSpan(), mapq);
    }

    std::fprintf(fpOut, "\ttp:A:%c", (ovl.IsSecondary ? 'S' : 'P'));

    std::fprintf(fpOut, "\tNS:i:%d", ovl.NumSeeds);
    std::fprintf(fpOut, "\tNM:i:%d\tIT:f:%.4lf\tSC:i:%d", ovl.EditDistance, 100.0 * identity,
                 static_cast<int32_t>(ovl.Score));
    std::fprintf(fpOut, "\tAT:Z:%s\tBT:Z:%s", AtypeStr.c_str(), BtypeStr.c_str());

    if (writeCigar) {
        std::fprintf(fpOut, "\tcg:Z:");
        if (ovl.Cigar.empty()) {
            std::fprintf(fpOut, "*");
        } else {
            for (const auto& op : ovl.Cigar) {
                std::fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
            }
        }
    }

    if (ovl.Avars.empty()) {
        std::fprintf(fpOut, "\tVQ:Z:*");
    } else {
        if (ovl.Arev) {
            auto vars = Pancake::ReverseComplement(ovl.Avars, 0, ovl.Avars.size());
            std::fprintf(fpOut, "\tVQ:Z:%s", vars.c_str());
        } else {
            std::fprintf(fpOut, "\tVQ:Z:%s", ovl.Avars.c_str());
        }
    }

    if (ovl.Bvars.empty()) {
        std::fprintf(fpOut, "\tVT:Z:*");
    } else {
        if (ovl.Brev) {
            auto vars = Pancake::ReverseComplement(ovl.Bvars, 0, ovl.Bvars.size());
            std::fprintf(fpOut, "\tVT:Z:%s", vars.c_str());
        } else {
            std::fprintf(fpOut, "\tVT:Z:%s", ovl.Bvars.c_str());
        }
    }

    std::fprintf(fpOut, "\n");
}

void OverlapWriterBase::PrintOverlapAsSAM(std::FILE* fpOut, const Overlap& ovl, const char* query,
                                          int64_t queryLen, const std::string& Aname,
                                          const std::string& Bname, bool writeIds, bool writeCigar)
{
    // double identity = static_cast<double>(ovl.Identity);
    // if (identity == 0.0 && ovl.EditDistance >= 0.0) {
    //     const double editDist = ovl.EditDistance;
    //     const double qSpan = ovl.ASpan();
    //     const double tSpan = ovl.BSpan();
    //     const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
    //     const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
    //     identity = std::min(identityQ, identityT);
    // }

    // The format specifies coordinates always in the FWD strand.
    const int32_t tStart = ovl.BstartFwd() + 1;
    const int32_t tIsRev = ovl.Brev;
    const std::string AtypeStr = OverlapTypeToStringSingleChar(ovl.Atype);
    const std::string BtypeStr = OverlapTypeToStringSingleChar(ovl.Btype);
    const int32_t flag = tIsRev ? 16 : 0;
    int32_t mapq = 60;
    const std::string seq =
        (ovl.Brev) ? ReverseComplement({query, static_cast<size_t>(queryLen)}, 0, queryLen)
                   : std::string(query, queryLen);
    const std::string qual = "*";

    // Query and target names, flag, pos and mapq.
    if (writeIds) {
        std::fprintf(fpOut, "%09d\t%d\t%09d\t%d\t%d", ovl.Aid, flag, ovl.Bid, tStart, mapq);
    } else {
        std::fprintf(fpOut, "%s\t%d\t%s\t%d\t%d", Aname.c_str(), flag, Bname.c_str(), tStart, mapq);
    }

    // Write the CIGAR only if specified, for speed.
    if (writeCigar) {
        if (ovl.Cigar.empty()) {
            std::fprintf(fpOut, "\t*");
        } else {
            if (ovl.Brev) {
                const std::string clipBack =
                    (ovl.Astart > 0) ? std::to_string(ovl.Astart) + "S" : "";
                const std::string clipFront =
                    (ovl.Aend < ovl.Alen) ? std::to_string(ovl.Alen - ovl.Aend) + "S" : "";
                std::fprintf(fpOut, "\t%s", clipFront.c_str());
                for (auto it = ovl.Cigar.rbegin(); it != ovl.Cigar.rend(); ++it) {
                    const auto& op = *it;
                    std::fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
                }
                std::fprintf(fpOut, "%s", clipBack.c_str());

            } else {
                const std::string clipFront =
                    (ovl.Astart > 0) ? std::to_string(ovl.Astart) + "S" : "";
                const std::string clipBack =
                    (ovl.Aend < ovl.Alen) ? std::to_string(ovl.Alen - ovl.Aend) + "S" : "";
                std::fprintf(fpOut, "\t%s", clipFront.c_str());
                for (const auto& op : ovl.Cigar) {
                    std::fprintf(fpOut, "%u%c", op.Length(), ConstexprTypeToChar(op.Type()));
                }
                std::fprintf(fpOut, "%s", clipBack.c_str());
            }
        }
    } else {
        std::fprintf(fpOut, "\t*");
    }

    std::fprintf(fpOut, "\t*\t0\t0\t%s\t%s", seq.c_str(), qual.c_str());
    std::fprintf(fpOut, "\tAT:Z:%s\tBT:Z:%s", AtypeStr.c_str(), BtypeStr.c_str());
    std::fprintf(fpOut, "\n");
}

std::string OverlapWriterBase::PrintOverlapAsM4(const Overlap& ovl, bool writeCigar)
{
    return PrintOverlapAsM4(ovl, "", "", true, writeCigar);
}

std::string OverlapWriterBase::PrintOverlapAsM4(const Overlap& ovl, const std::string& Aname,
                                                const std::string& Bname, bool writeIds,
                                                bool writeCigar)
{
    double identity = static_cast<double>(ovl.Identity);
    if (identity == 0.0 && ovl.EditDistance >= 0.0) {
        const double editDist = ovl.EditDistance;
        const double qSpan = ovl.ASpan();
        const double tSpan = ovl.BSpan();
        const double identityQ = (qSpan != 0) ? ((qSpan - editDist) / qSpan) : -2.0;
        const double identityT = (tSpan != 0) ? ((tSpan - editDist) / tSpan) : -2.0;
        identity = std::min(identityQ, identityT);
    }

    // The format specifies coordinates always in the FWD strand.
    const int32_t tStart = ovl.BstartFwd();
    const int32_t tEnd = ovl.BendFwd();
    const int32_t tIsRev = ovl.Brev;
    const int32_t tLen = ovl.Blen;
    const std::string AtypeStr = OverlapTypeToString(ovl.Atype);

    char buffA[100], buffB[100];
    char idtBuff[100];
    std::sprintf(buffA, "%09d", ovl.Aid);
    std::sprintf(buffB, "%09d", ovl.Bid);
    std::sprintf(idtBuff, "%.2lf", 100.0 * identity);

    std::ostringstream oss;
    if (writeIds) {
        oss << buffA << " " << buffB;
    } else {
        oss << Aname << " " << Bname;
    }

    std::string cigar;
    std::string cigarSep;
    if (writeCigar) {
        cigarSep = " ";
        cigar = (ovl.Cigar.empty()) ? "*" : ovl.Cigar.ToStdString();
    }

    oss << " " << static_cast<int32_t>(ovl.Score) << " " << idtBuff << " "
        << static_cast<int32_t>(ovl.Arev) << " " << ovl.Astart << " " << ovl.Aend << " " << ovl.Alen
        << " " << static_cast<int32_t>(tIsRev) << " " << tStart << " " << tEnd << " " << tLen << " "
        << AtypeStr << cigarSep << cigar;

    return oss.str();
}

std::ostream& operator<<(std::ostream& os, const Overlap& ovl)
{
    os << OverlapWriterBase::PrintOverlapAsM4(ovl, "", "", true, false);
    return os;
}

}  // namespace Pancake
}  // namespace PacBio
