#include <iostream>
#include <fstream>
#include "FileUtil.hpp"
#include "common.hpp"
using namespace std;

bool getECOrder(string ec, vector<string> &order,
        unordered_set<string> &ecSet) {
   ifstream in(ec);
   if (!in.is_open()) { return false; }
   string inp;
   while (getline(in, inp)) {
        string ec = parseString(inp, "\t", 0)[1];
        order.push_back(ec);
        ecSet.emplace(ec);
   }
   return true;
}

bool readTranscriptomeFile(string file, int &transcript_count,
        unordered_map<string, int> &indexMap) {
    ifstream in(file);
    if (!in.is_open()) { return false; }

    string inp;
    while(getline(in, inp)) {
        if (inp.size() == 0 || inp[0] != '>') { continue; }
        int start = inp.find(TRANSCRIPTOME_START);
        if (start == string::npos) {
            ++transcript_count;
            continue;
        }
        start += string(TRANSCRIPTOME_START).size();
        int end = inp.find(TRANSCRIPTOME_END, start);
        if (end == string::npos) {
            end = inp.size();
        }
        indexMap.emplace(inp.substr(start, end - start), transcript_count);
        ++transcript_count;
    }

    in.close();
    return true;
}

bool readTranscriptome(vector<string> &files,
        unordered_map<string, int> &indexMap) {
    int transcript_count = 0;
    for (int i = 0; i < files.size(); ++i) {
        if (!readTranscriptomeFile(files[i], transcript_count, indexMap)) {
            cerr << "  WARNING: unable to read " << files[i] << endl;
        }
    }
    return true;
}

