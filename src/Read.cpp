#include <algorithm> /* sort, set_intersection, unique */
#include "Read.hpp"
using namespace std;

Read::Read() {}

Read::Read(const seqan::BamAlignmentRecord &alignment) {
    seen[0] = seen[1] = 0;
    NH[0] = NH[1] = (seqan::hasFlagMultiple(alignment)) ? -1 : 0;
    fillNH(alignment);
    for (int i = 0; i < 4; ++i) {
        ECs[i] = new vector<int>;
    }
}

Read::~Read() {
    for (int i = 0; i < 4; ++i) {
        delete ECs[i];
    }
}

int Read::needsNH() {
    if (NH[0] == -1) { return 0; }
    else if (NH[1] == -1) { return 1; }
    else { return -1; }
}

void Read::fillNH(const seqan::BamAlignmentRecord &alignment) {
    seqan::BamTagsDict tags(alignment.tags);
    int nh = 0, id;
    seqan::findTagKey(id, tags, "NH") && seqan::extractTagValue(nh, tags, id);
    if (seqan::hasFlagLast(alignment)) {
        NH[1] = nh;
    } else {
        NH[0] = nh;
    }
}

void Read::addEC(vector<int> &EC, bool first, bool reverse) {
    int i;
    if (first) {
        ++seen[0];
        i = 0;
    } else {
        ++seen[1];
        i = 1;
    }
    if (reverse) { i += 2; }
    ECs[i]->insert(ECs[i]->end(), EC.begin(), EC.end());
}

bool Read::isComplete() {
    return NH[0] == seen[0] && NH[1] == seen[1];
}

string Read::getEC(bool paired, bool genomebam) {
    vector<int> EC;
    for (int i = 0; i < 4; ++i) {
        sort(ECs[i]->begin(), ECs[i]->end());
    }
    if (paired && (!genomebam || (ECs[0]->size() + ECs[2]->size() != 0
                            && ECs[1]->size() + ECs[3]->size() != 0))) {
        set_intersection(ECs[0]->begin(), ECs[0]->end(),
            ECs[3]->begin(), ECs[3]->end(),
            back_inserter(EC));
        set_intersection(ECs[1]->begin(), ECs[1]->end(),
            ECs[2]->begin(), ECs[2]->end(),
            back_inserter(EC));
    } else {
        if (genomebam && ECs[0]->size() + ECs[2]->size() == 0) {
            vector<int> *temp = ECs[0];
            ECs[0] = ECs[1];
            ECs[1] = temp;
            temp = ECs[2];
            ECs[2] = ECs[3];
            ECs[3] = temp;
        }
        for (int i = 0; i < 3; i += 2) {
            EC.insert(EC.end(), ECs[i]->begin(), ECs[i]->end());
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
