#ifndef __FILE_META_INFO__
#define __FILE_META_INFO__

struct FileMetaInfo {
    int fileNum, start, end, count;
    FileMetaInfo(int fileNum, int start, int end, int count) :
        fileNum(fileNum), start(start), end(end), count(count) {};
};

#endif
