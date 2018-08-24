#include "SequenceNumID.hpp"

int getTranscriptID(SequenceNumID n) {
    return (int) (n / SHIFT);
}

int getExonIndex(SequenceNumID n) {
    return (int) (n % SHIFT);
}

SequenceNumID getSequenceNumID(int transcript, int exon) {
    if (exon >= SHIFT || transcript >= MAX_TRANSCRIPT) {    
        return -1;
    }
    return transcript * SHIFT + exon;
}

bool compareSequenceNumID(const SequenceNumID &a, const SequenceNumID &b) {
    return getTranscriptID(a) < getTranscriptID(b);
}

