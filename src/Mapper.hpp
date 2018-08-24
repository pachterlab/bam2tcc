#ifndef __MAPPER_HPP__
#define __MAPPER_HPP__

#include <vector>
#include <unordered_map>
#include <deque>
#include <set>
#include <mutex>
#include <condition_variable>
#include "TCC_Matrix.hpp"
#include "ChromMetaInfo.hpp"
#include "Read.hpp"
#include "Transcript.hpp"
#include "Semaphore.hpp"

#define TRANSCRIPT_ID_TAG "transcript_id"

class Mapper {
private:
    std::vector<std::string> gffs;
    std::vector<std::string> sams;
    std::unordered_map<std::string, int> *indexMap;
    std::unordered_map<std::string, ChromMetaInfo> *chroms;
    std::unordered_map<std::string, Read*> *reads;
    Semaphore *readsSem;
    TCC_Matrix *matrix;
    bool paired;
    
    bool readGFF(ChromMetaInfo &inf, std::deque<Transcript> &chrom);
    bool readSAM(ChromMetaInfo &inf, std::deque<Transcript> &chrom,
            bool genomebam, bool rapmap, bool sameQName);
    bool mapToChrom(ChromMetaInfo &inf,
            bool genomebam, bool rapmap, bool sameQName,
            int thread, std::condition_variable &cv, std::mutex &m,
            std::queue<int> &completed);
    bool getChromsGFF(int filenumber, int &transcriptCount);
    bool getChromsSAM(int filenumber);
    bool getChromsGFFs();
    bool getSameQName(int filenumber, bool &same);
    std::string getSamPGName(int filenumber);
    bool getPG(int filenumber, bool &genomebam, bool &rapmap);
    bool writeCellsFiles(std::string outprefix);
public:
    Mapper(std::vector<std::string> gffs, std::vector<std::string> sams,
            std::vector<std::string> fas, bool paired);
    ~Mapper();
    bool mapReads(int nThreads);
    bool writeToFile(std::string outprefix, bool full, std::string ec);
};
#endif

