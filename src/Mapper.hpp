#ifndef __MAPPER_HPP__
#define __MAPPER_HPP__

#include <vector>
#include <unordered_map>
#include <deque>
#include <set>
#include <mutex>
#include <condition_variable>
#include "TCC_Matrix.hpp"
#include "FileMetaInfo.hpp"
#include "Read.hpp"
#include "Transcript.hpp"
#include "Semaphore.hpp"

#define DEBUG 0
#define READ_DIST 0
#define TRANSCRIPT_ID_TAG "transcript_id"

class Mapper {
private:
    std::vector<std::string> gffs;
    std::vector<std::string> sams;
    std::unordered_map<std::string, int> *indexMap;
    std::unordered_map<std::string, FileMetaInfo> chroms;
    std::vector<std::unordered_map<std::string, Read*>*> reads;
    std::vector<Semaphore*> readsSems;
    std::vector<std::unordered_set<std::string>*> unmappedQNames;
    std::vector<Semaphore*> unmappedQNamesSems;
    TCC_Matrix *matrix;
    bool paired, recordUnmapped, pgProvided, genomebam, rapmap;
#if READ_DIST
    std::vector<std::unordered_set<std::string>*> mappedQNames;
    std::vector<Semaphore*> mappedQNamesSems;
#endif
#if DEBUG
    Semaphore debugOutSem;
#endif
    
    bool readGFF(FileMetaInfo &inf, std::deque<Transcript> &chrom);
    bool readSAM(FileMetaInfo &inf, std::deque<Transcript> &chrom,
            bool genomebam, bool rapmap, bool sameQName);
    bool mapToChrom(FileMetaInfo &gffInf, FileMetaInfo samInf,
            bool genomebam, bool rapmap, bool sameQName,
            int thread, std::condition_variable &cv, std::mutex &m,
            std::queue<int> &completed);
    bool getChromsGFF(int filenumber, int &transcriptCount);
    bool getChromsSAM(int filenumber,
            std::unordered_map<std::string, FileMetaInfo> &inf);
    bool getChromsGFFs();
    bool getSameQName(int filenumber, bool &same);
    std::string getSamPGName(int filenumber);
    bool getPG(int filenumber, bool &genomebam, bool &rapmap);
    bool mapUnmapped(int samNum, int start, int end, bool genomebam);
    bool writeCellsFiles(std::string outprefix);
    bool writeUnmapped(std::vector<std::string> &unmappedOut);
#if READ_DIST
    bool writeMapped(std::vector<std::string> &mappedOut);
#endif
public:
    Mapper(std::vector<std::string> gffs, std::vector<std::string> sams,
            std::vector<std::string> fas, bool paired, bool recordUnmapped,
            bool pgProvided, bool genomebam, bool rapmap);
    ~Mapper();
    bool mapReads(int nThreads);
    bool writeToFile(std::string outprefix,
            std::vector<std::string> &unmappedOut,
#if READ_DIST
            std::vector<std::string> &mappedOut,
#endif
            bool full, std::string ec);
};
#endif

