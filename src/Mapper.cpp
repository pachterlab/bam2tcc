#include <seqan/gff_io.h>
#include <seqan/bam_io.h>
#include <fstream>
#include <future>
#include "Mapper.hpp"
#include "Exon.hpp"
#include "FileUtil.hpp"
#include "common.hpp"
using namespace std;

Mapper::Mapper(vector<string> gffs, vector<string> sams, vector<string> fas,
        bool paired) : gffs(gffs), sams(sams), paired(paired) {
    indexMap = new unordered_map<string, int>;
    chroms = new unordered_map<string, ChromMetaInfo>;
    reads = new unordered_map<string, Read*>;
    readsSem = new Semaphore();
    matrix = new TCC_Matrix(sams.size());
#if DEBUG
    debugOutSem = new Semaphore();
#endif

    readTranscriptome(fas, *indexMap);
    getChromsGFFs();
}

Mapper::~Mapper() {
    delete indexMap;
    delete chroms;
    for (auto it = reads->begin(); it != reads->end(); ++it) {
        delete it->second;
    }
    delete reads;
    delete readsSem;
    delete matrix;
}

bool Mapper::readGFF(ChromMetaInfo &inf, deque<Transcript> &chrom) {
    seqan::GffFileIn gff;
    if (!seqan::open(gff, gffs[inf.gffNum].c_str())) { return false; }
    int line = 0, transcriptCount = inf.transcriptCount;
    seqan::GffRecord rec;
    while (line < inf.gffStart) {
        ++line;
        seqan::readRecord(rec, gff);
    }

    multiset<Transcript> transcripts;
    Transcript transcript;
    string tofind = TRANSCRIPT_ID_TAG;

    while (true) {
        string type = lower(seqan::toCString(rec.type));
        if (type.compare("transcript") == 0) {
            transcripts.insert(transcript);

            int id = transcriptCount;
            if (indexMap->size()) {
                string transcript_id;
                for (int i = 0; i < length(rec.tagNames); ++i) {
                    if (tofind.compare(seqan::toCString(rec.tagNames[i])) == 0)
                    {
                        transcript_id = seqan::toCString(rec.tagValues[i]);
                        break;
                    }
                }
                if (indexMap->find(transcript_id) == indexMap->end()) {
                    id = -1;
                } else {
                    id = indexMap->at(transcript_id);
                }
            }
            transcript = Transcript(id, rec);
            ++transcriptCount;
        } else if (type.compare("exon") == 0) {
            transcript.addExonEntry(rec);
        }
        ++line;
        if (line == inf.gffEnd) { break; }
        seqan::readRecord(rec, gff);
    }

    for (auto it = transcripts.begin(); it != transcripts.end(); ++it) {
        chrom.push_back(*it);
    }

    return true;
}

vector<Exon> getAlignmentExons(const seqan::BamAlignmentRecord &alignment) {
    vector<Exon> exons;
    int start = alignment.beginPos, end = start;
    for (uint i = 0; i < seqan::length(alignment.cigar); ++i) {
        switch (alignment.cigar[i].operation) {
            case 'M':
            case 'D':
            case '=':
            case 'X': end += alignment.cigar[i].count;
                      break;
            case 'N': exons.push_back(Exon(start, end));
                      start = end + alignment.cigar[i].count;
                      end = start;
                      break;
            default: /* do nothing */ break;
        }
    }
    exons.push_back(Exon(start, end));
    return exons;
}

bool Mapper::readSAM(ChromMetaInfo &inf, deque<Transcript> &chrom,
        bool genomebam, bool rapmap, bool sameQName) {
    seqan::BamFileIn bam;
    if (!seqan::open(bam, sams[inf.samNum].c_str())) { return false; }
    int line = 0;
    seqan::BamAlignmentRecord rec;
    seqan::BamHeader head;
    seqan::readHeader(head, bam);
    while (line < inf.samStart) {
        ++line;
        readRecord(rec, bam);
    }

    while (true) {
        vector<int> EC;
        if (rapmap) {
            EC = {rec.rID};
        } else {
            while (!chrom.empty()
                && chrom.front().getEnd() <= rec.beginPos) {
                chrom.pop_front();
            }
            if (chrom.empty()) { return true; }

            if (!seqan::hasFlagUnmapped(rec)
                    && ((!genomebam && (!seqan::hasFlagMultiple(rec)
                        || (seqan::hasFlagAllProper(rec)
                            && rec.rID == rec.rNextId)))
                    || (genomebam && (!paired || (rec.rID == rec.rNextId
                            && seqan::hasFlagMultiple(rec))))))
            {
                vector<Exon> alignmentExons = getAlignmentExons(rec);
                for (auto it = chrom.begin(); it != chrom.end(); ++it) {
                    if (it->mapsToTranscript(alignmentExons, genomebam)) {
                        EC.push_back(it->getID());
                    }
                }
            }
        }

        string qName = seqan::toCString(rec.qName);
        if (!sameQName) {
            qName = qName.substr(0, qName.size() - 2);
        }

        readsSem->dec();
        Read *read;
        if (reads->find(qName) == reads->end()) {
            read = new Read(rec);
            reads->emplace(qName, read);
        } else {
            read = reads->at(qName);
        }
        int needsNH = read->needsNH();
        if (needsNH != -1 && (bool) needsNH == seqan::hasFlagLast(rec)) {
            read->fillNH(rec);
        }
        read->addEC(EC,
                (!seqan::hasFlagMultiple(rec) || seqan::hasFlagFirst(rec)),
                seqan::hasFlagRC(rec));
        bool complete = read->isComplete();
        if (complete) {
            reads->erase(qName);
        }
        readsSem->inc();

        if (complete) {
            string stringEC = read->getEC(paired, genomebam);
            if (stringEC.size() == 0) { /* Do something */ }
            else {
                matrix->inc_TCC(stringEC, inf.samNum);
            }
            delete read;
        }

        ++line;
#if DEBUG
        cout << "." << flush;
#endif
        if (line == inf.samEnd) { break; }
        readRecord(rec, bam);
    }

    return true;
}

bool Mapper::mapToChrom(ChromMetaInfo &inf,
        bool genomebam, bool rapmap, bool sameQName,
        int thread, condition_variable &cv, mutex &m, queue<int> &completed) {
    deque<Transcript> *chrom = new deque<Transcript>;
#if DEBUG
    debugOutSem->dec();
    cout << "    Thread " << thread << " reading GFF" << endl;
    debugOutSem->inc();
#endif
    if (!rapmap && !readGFF(inf, *chrom)) { return false; }
#if DEBUG
    debugOutSem->dec();
    cout << "    Thread " << thread << " reading SAM" << endl;
    debugOutSem->inc();
#endif
    if (!readSAM(inf, *chrom, genomebam, rapmap, sameQName)) { return false; }
    delete chrom;

    m.lock();
    completed.push(thread);
    m.unlock();
    cv.notify_one();

#if DEBUG
    debugOutSem->dec();
    cout << "    Thread " << thread << " complete" << endl;
    debugOutSem->inc();
#endif

    return true;
}       

bool Mapper::getChromsGFF(int filenumber, int &transcriptCount) {
    ifstream in(gffs[filenumber]);
    if (!in.is_open()) { return false; }
    string inp, currChrom;
    int line = 0, start;
    while(getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '#') { continue; }
        ++line;
        vector<string> parsed = parseString(inp, "\t", 3);
        string chrom = parsed[0], type = parsed[2];
        if (lower(type).compare("transcript") != 0) { continue; }
        currChrom = chrom; 
        start = line;
        break;
    } 
    while(getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '#') { continue; }
        ++line;
        vector<string> parsed = parseString(inp, "\t", 3);
        string chrom = parsed[0], type = parsed[2];
        if (lower(type).compare("transcript") != 0) { ++transcriptCount; }
        if (chrom.compare(currChrom) != 0) {
            chroms->emplace(currChrom, ChromMetaInfo(filenumber, start, line,
                        transcriptCount));
            start = line;
            currChrom = chrom;
        }
    }
    chroms->emplace(currChrom, ChromMetaInfo(filenumber, start, line + 1,
                transcriptCount));
    ++transcriptCount;
    in.close();
    return true; 
}

bool Mapper::getChromsSAM(int filenumber) {
    ifstream in(sams[filenumber]);
    if (!in.is_open()) { return false; }

    for (auto it = chroms->begin(); it != chroms->end(); ++it) {
        it->second.clearSAMInfo();
    }

    string inp, currChrom;
    unordered_set<string> notInGFF;
    int line = 0, start;
    while(getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '@') { continue; }
        ++line;
        string chrom = parseString(inp, "\t", 3)[2];
        if (chrom.compare("*") == 0) { continue; }
        currChrom = chrom;
        start = line;
        break;
    } 
    while(getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '@') { continue; }
        ++line;
        string chrom = parseString(inp, "\t", 3)[2];
        if (chrom.compare("*") == 0) { continue; }
        if (chrom.compare(currChrom) != 0) {
            auto it = chroms->find(currChrom);
            if (it == chroms->end()) {
                notInGFF.emplace(currChrom);
            } else {
                it->second.setSAM(filenumber, start, line);
            }
            start = line;
            currChrom = chrom;
        }
    }
    auto it = chroms->find(currChrom);
    if (it == chroms->end()) {
        notInGFF.emplace(currChrom);
    } else {
        it->second.setSAM(filenumber, start, line + 1);
    }

    if (!notInGFF.empty()) {
        cerr << "  WARNING: " << notInGFF.size() << " chromosomes/scaffolds "
            << "present in " << sams[filenumber] << " but not in GFFs:" << endl;
        auto it = notInGFF.begin();
        cerr << *it << flush;
        while (it != notInGFF.end()) {
            cerr << ',' << *it << flush;
            ++it;
        }
        cerr << endl;
    }

    in.close();
    return true; 
}

bool Mapper::getChromsGFFs() {
    int transcript_count = 0;
    for (int i = 0; i < gffs.size(); ++i) {
        if (!getChromsGFF(i, transcript_count)) {
            cerr << "WARNING: error while reading " << gffs[i] << endl;
        }
    }
    return true;
}

bool Mapper::getSameQName(int filenumber, bool &same) {
    ifstream in(sams[filenumber]);
    if (!in.is_open()) { return false; }
    same = true;
    bool one_seen = false, two_seen = false;
    string inp;
    while (getline(in, inp)) {
        string qName = parseString(inp, "\t", 1)[0];
        if (qName.size() < 2) {
            return true;;
        }
        if (!isdigit(qName[qName.size() - 2])) {
            if (qName[qName.size() - 1] == '1') {
                if (one_seen && two_seen) {
                    same = false;
                    break;
                } else {
                    one_seen = true;
                }
            } else if (qName[qName.size() - 1] == '2') {
                two_seen = true;
            } else {
                break;
            }
        } else {
            break;
        }
    }
    in.close();
    return true;
}

string Mapper::getSamPGName(int filenumber) {
    ifstream in(sams[filenumber]);
    if (!in.is_open()) {
        return "";
    }
    string inp, pg = "N/A";
    while (getline(in, inp)) {
        if (inp.size() != 0 && inp.substr(0, 3).compare("@PG") == 0) {
            int start = inp.find("ID:");
            if (start == string::npos) { break; }
            start += 3;
            int end = inp.find('\t', start);
            if (end == string::npos) { break; }
            pg = inp.substr(start, end - start);
            break;
        } else if (inp.size() != 0 && inp[0] != '@') {
            break;
        }
    }
    in.close();
    return pg;
}

bool Mapper::getPG(int filenumber, bool &genomebam, bool &rapmap) {
    string pg = getSamPGName(filenumber);
    if (pg.size() == 0) { return false; }
    genomebam = false, rapmap = false;
    if (pg.compare("kallisto") == 0) { genomebam = true; }
    else if (pg.compare("rapmap") == 0) { rapmap = true; }
    return true;
}

bool Mapper::mapReads(int nThreads) {
    if (nThreads <= 0) {
        cerr << "  ERROR: cannot run with nonpositive number of threads."
            << endl;
        return false;
    }

    bool genomebam = false, rapmap = false, sameQName = false;

    condition_variable cv;
    mutex m;
    queue<int> completed;
    future<bool> threads[nThreads];

    for (int i = 0; i < nThreads; ++i) {
        completed.push(i);
    }

    for (int i = 0; i < sams.size(); ++i) {
#if DEBUG
        debugOutSem->dec();
        cout << "  Mapping " << sams[i] << endl;
        debugOutSem->inc();
#endif
        if (!getChromsSAM(i)) {
            cerr << "  WARNING: error while reading " << sams[i] << endl;
            continue;
        }

        for (auto chrom = chroms->begin(); chrom != chroms->end(); ++chrom) {
            if (!chrom->second.isSAMSet()) { continue; }
            if (getSameQName(chrom->second.samNum, sameQName)
                    && getPG(chrom->second.samNum, genomebam, rapmap)) {
                unique_lock<mutex> lk(m);
                if (completed.empty()) {
                    cv.wait(lk, [&completed] { return !completed.empty(); });
                }
                int done = completed.front();
                completed.pop();
                lk.unlock();
                if (threads[done].valid()) {
                    if (!threads[done].get()) {
                        cerr << "  WARNING: thread failed." << endl;
                    }
                }
                threads[done] = async(launch::async, &Mapper::mapToChrom, this,
                        ref(chrom->second), genomebam, rapmap, sameQName,
                        done, ref(cv), ref(m), ref(completed));
#if DEBUG
                debugOutSem->dec();
                cout << "    Thread " << done << " launched on "
                   <<  chrom->first << endl;
                debugOutSem->inc();
#endif
            } else {
                cerr << "  WARNING: error reading chromosome "
                    << chrom->first << endl;
            }
        }
    }

    for (int i = 0; i < nThreads; ++i) {
        if (threads[i].valid()) {
            if (!threads[i].get()) {
                    cerr << "  WARNING: thread failed." << endl;
            }
        }
    }

    cerr << reads->size() << " reads not completed." << flush;
    if (genomebam && sams.size() == 1) {
        cerr << " Placing in TCC." << flush;
        for (auto it = reads->begin(); it != reads->end(); ++it) {
            string stringEC = it->second->getEC(paired, genomebam);
            if (stringEC.size() == 0) { /* Do something? */ }
            else {
                matrix->inc_TCC(stringEC, 0);
            }
            delete it->second;
        }
        reads->clear();
    }
    cerr << endl;

    return true;
}

bool Mapper::writeCellsFiles(string outprefix) {
    ofstream out(outprefix + ".cells");
    if (!out.is_open()) { return false; }
    for (auto file = sams.begin(); file != sams.end(); ++file) {
        int start = file->find_last_of('/');
        if (start == string::npos) { start = -1; }
        ++start;

        /* TODO: stop being a cancer upon this world. */
        int end = file->find_last_of(".bam");
        if (end == string::npos) {
            end = file->find_last_of(".sam");
            if (end == string::npos) { end = file->size() + 3; }
        }
        end -= 3;

        if (end > start) { out << file->substr(start, end - start) << endl; }
        else { out << *file << endl; }
    }
    out.close();
    return true;
}

bool Mapper::writeToFile(string outprefix, bool full, string ec) {
    if (ec.size() == 0) {
        if (full) { matrix->write_to_file(outprefix); }
        else { matrix->write_to_file_sparse(outprefix); }
    } else {
        vector<string> order;
        unordered_set<string> ecSet;
        getECOrder(ec, order, ecSet);
        if (full) { matrix->write_to_file_in_order(outprefix, order, ecSet); }
        else { matrix->write_to_file_in_order_sparse(outprefix, order, ecSet); }
    }
    writeCellsFiles(outprefix);
    return true;
}

