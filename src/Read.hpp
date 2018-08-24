#ifndef __READ_HPP__
#define __READ_HPP__

#include <vector>
#include <seqan/bam_io.h>

class Read {
private:
    int NH[2];
    int seen[2];
    std::vector<int> *ECs[4];
public:
    Read();
    Read(const seqan::BamAlignmentRecord &alignment);
    ~Read();
    int needsNH();
    void fillNH(const seqan::BamAlignmentRecord &alignment);
    void addEC(std::vector<int> &EC, bool first, bool reverse);
    bool isComplete();
    std::string getEC(bool paired, bool genomebam=false);
};

#endif
