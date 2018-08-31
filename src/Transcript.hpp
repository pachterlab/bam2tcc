#ifndef __TRANSCRIPT_HPP__
#define __TRANSCRIPT_HPP__

#include <seqan/gff_io.h>
#include "Exon.hpp"

class Transcript {
private:
    int id, start, end;
    std::vector<Exon> exons;
public:
    Transcript();
    Transcript(int id, const seqan::GffRecord &entry);
    int getID();
    int getEnd();
    void addExonEntry(const seqan::GffRecord &entry);
    bool mapsToTranscript(const std::vector<Exon> &alignmentExons,
            bool genomebam);
    bool operator<(const Transcript &other) const;
};

#endif
