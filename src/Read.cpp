#include <algorithm> /* sort, set_intersection, unique */
#include "Read.hpp"
using namespace std;

Read::Alignment::Alignment(int rName, int rNext, int pos, int nextPos,
        bool first, bool reverse, const vector<int> &EC) : rName(rName),
    rNext(rNext), pos(pos), nextPos(nextPos), first(first), reverse(reverse),
    EC(EC) {}

Read::Alignment::~Alignment() {}

Read::Pair::Pair(const vector<int> &EC) : EC1(EC) {}

Read::Pair::Pair(const vector<int> &EC1, const vector<int> &EC2) :
    EC1(EC1), EC2(EC2) {}

Read::Pair::~Pair() {}

Read::Read() {}

Read::Read(const seqan::BamAlignmentRecord &alignment, const vector<int> &EC) {
    paired = true;
    seen[0] = 0;
    seen[1] = 0;
    NH[0] = -1;
    NH[1] = -1;

    if (!seqan::hasFlagMultiple(alignment)) {
        paired = false;
        NH[1] = 0;
    }
    addAlignment(alignment, EC);
}

Read::~Read() {}

int Read::getNH(const seqan::BamAlignmentRecord &alignment) {
    if (seqan::hasFlagUnmapped(alignment)) {
        return 1;
    }
    seqan::BamTagsDict tags(alignment.tags);
    int nh = 1, id;
    seqan::findTagKey(id, tags, "NH") && seqan::extractTagValue(nh, tags, id);
    return nh;
}

void Read::addAlignment(const seqan::BamAlignmentRecord &alignment,
           const vector<int> &EC) {
    int i = (!paired || seqan::hasFlagFirst(alignment)) ? 0 : 1;
    ++seen[i];
    if (NH[i] == -1) {
        NH[i] = getNH(alignment);
    }
   
    if (!paired) {
        pairs.emplace_back(EC);
        return;
    }

    auto a2 = alignments.begin();
    while (a2 != alignments.end()) {
        if (alignment.rID == a2->rName && alignment.rNextId == a2->rNext
                && alignment.beginPos == a2->nextPos
                && alignment.pNext == a2->pos
                && seqan::hasFlagFirst(alignment) != a2->first) {
            break;
        }
        ++a2;
    }

    if (a2 == alignments.end()) {
        alignments.emplace_back(Alignment(alignment.rID, alignment.rNextId,
                    alignment.beginPos, alignment.pNext,
                    seqan::hasFlagFirst(alignment),
                    seqan::hasFlagRC(alignment), EC));
    } else {
        if (seqan::hasFlagRC(alignment) != a2->reverse) {
            pairs.emplace_back(a2->EC, EC);
        }
        alignments.erase(a2);
    }
}

bool Read::isComplete() {
    return NH[0] == seen[0] && NH[1] == seen[1];
}

string Read::getEC(bool genomebam) {
    vector<int> EC;
    for (auto p = pairs.begin(); p != pairs.end(); ++p) {
        if (paired && (!genomebam
                        || p->EC1.size() + p->EC2.size() != 0)) {
            sort(p->EC1.begin(), p->EC1.end());
            sort(p->EC2.begin(), p->EC2.end());
            set_intersection(p->EC1.begin(), p->EC1.end(),
                    p->EC2.begin(), p->EC2.end(),
                    back_inserter(EC));
        } else {
            if (genomebam && p->EC1.size() == 0) {
                p->EC1 = p->EC2;
            }
            EC.insert(EC.end(), p->EC1.begin(), p->EC1.end());
        }
    }
    sort(EC.begin(), EC.end());
    EC.erase(unique(EC.begin(), EC.end()), EC.end());
    if (EC.size() == 0) { return ""; }
    string stringEC = to_string(EC[0]);
    for (int i = 1; i < EC.size(); ++i) {
        stringEC += ',' + to_string(EC[i]);
    }
    return stringEC;
}
