#ifndef __FILE_UTIL_HPP__
#define __FILE_UTIL_HPP__

#include <vector>
#include <string>
#include <unordered_map>
#include <unordered_set>

#define TRANSCRIPTOME_START ">"
#define TRANSCRIPTOME_END " "

bool readTranscriptome(std::vector<std::string> &files,
        std::unordered_map<std::string, int> &indexMap);

bool getECOrder(std::string ec, std::vector<std::string> &order,
            std::unordered_set<std::string> &ecSet);

#endif
