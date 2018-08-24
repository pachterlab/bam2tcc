#ifndef __CHROM_META_INFO__
#define __CHROM_META_INFO__

struct ChromMetaInfo {
    int gffNum, gffStart, gffEnd, transcriptCount;
    int samNum, samStart, samEnd;
    ChromMetaInfo(int gffNum, int gffStart, int gffEnd, int transcriptCount) :
        gffNum(gffNum), gffStart(gffStart), gffEnd(gffEnd),
        transcriptCount(transcriptCount), samNum(-1) {};
    void setSAM(int samNum, int samStart, int samEnd) {
        this->samNum = samNum;
        this->samStart = samStart;
        this->samEnd = samEnd;
    }
    bool isSAMSet() { return samNum != -1; }
    void clearSAMInfo() { samNum = -1; }
};

#endif
