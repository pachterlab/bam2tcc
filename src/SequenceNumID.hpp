#ifndef __SEQUENCE_NUM_ID_HPP__
#define __SEQUENCE_NUM_ID_HPP__

#include <cstdint>

#define SHIFT 10000
#define MAX_TRANSCRIPT 429496 /* ...I think. Make this less cancer later. */

typedef uint32_t SequenceNumID;
int getTranscriptID(SequenceNumID n);
int getExonIndex(SequenceNumID n);
SequenceNumID getSequenceNumID(int transcript, int exon);
bool compareSequenceNumID(const SequenceNumID &a, const SequenceNumID &b);

#endif

