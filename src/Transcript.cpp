#include "Transcript.hpp"
using namespace std;

Transcript::Transcript() {}

Transcript::Transcript(int id, const seqan::GffRecord &entry) : id(id),
        start(entry.beginPos),
        end(entry.endPos) {};

int Transcript::getID() { return id; }

int Transcript::getEnd() { return end; }

void Transcript::addExonEntry(const seqan::GffRecord &entry) {
    if (entry.strand == '+') {
        exons.push_back(Exon(entry.beginPos, entry.endPos));
    } else if (entry.strand == '-') {
        exons.insert(exons.begin(), Exon(entry.beginPos, entry.endPos));
    }
}

bool Transcript::mapsToTranscript(const vector<Exon> &alignmentExons,
        bool genomebam) {
    if (alignmentExons.begin()->start < start
            || alignmentExons[alignmentExons.size() - 1].end > end
            || alignmentExons.size() > exons.size()) {
        return false;
    }

    auto alignmentExon = alignmentExons.begin();
    bool aligning = false;
    for (auto exon = exons.begin(); exon != exons.end(); ++exon) {
        if (alignmentExon == alignmentExons.end()) { return true; }
        if (exon->start <= alignmentExon->start
                && alignmentExon->end <= exon->end
                && (genomebam
                    || ((alignmentExon == alignmentExons.begin()
                        || exon->start == alignmentExon->start)
                    && (alignmentExon == alignmentExons.end() - 1
                        || exon->end == alignmentExon->end)))) {
            aligning = true;
            ++alignmentExon;
        } else if (aligning) {
            return false;
        }
    }
    return alignmentExon == alignmentExons.end();
}

bool Transcript::operator<(const Transcript &other) const {
    return end < other.end;
}

