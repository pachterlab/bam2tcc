#ifndef __READ_HPP__
#define __READ_HPP__

#include <vector>
#include <seqan/bam_io.h>

class Read {
private:
    struct Alignment {
        int rName;
        int rNext;
        int pos;
        int nextPos;
        bool first;
        bool reverse;
        std::vector<int> EC;
        Alignment(int rName, int rNext, int pos, int nextPos, bool first,
                bool reverse, const std::vector<int> &EC);
        ~Alignment();
    };
    struct Pair {
        std::vector<int> EC1;
        std::vector<int> EC2;
        Pair(const std::vector<int> &EC);
        Pair(const std::vector<int> &EC1, const std::vector<int> &EC2);
        ~Pair();
    };
    bool paired;
    int NH[2];
    int seen[2];
    std::vector<Alignment> alignments;
    std::vector<Pair> pairs;
    int getNH(const seqan::BamAlignmentRecord &alignment);
public:
    Read();
    Read(const seqan::BamAlignmentRecord &alignment,
            const std::vector<int> &EC);
    ~Read();
    void addAlignment(const seqan::BamAlignmentRecord &alignment,
            const std::vector<int> &EC, bool genomebam);
    bool isComplete();
    std::string getEC(bool genomebam=false);
};

#endif
