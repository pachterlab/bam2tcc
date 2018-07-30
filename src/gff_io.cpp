#include <seqan/gff_io.h>
#include "gff_io.hpp"
#include "util.hpp"
#include "kallisto_util.hpp"
using namespace std;

/**
 * @brief Returns a vector allocated on the heap containing all exons in map.
 *
 * @param m     Map containing exons to place in vector.
 * @return      Vector containing exons.
 */
vector<Exon> *map_values(unordered_map<string, Exon> &m) {
    vector<Exon> *vec = new vector<Exon>;
    for (unordered_map<string, Exon>::iterator it = m.begin();
         it != m.end(); ++it) {
        vec->push_back(Exon(it->second));
    }
    return vec;
}

/**
 * @brief Populates and returns a newly allocated Sequence struct based on
 * information provided in string info, a line from a GTF file.
 *
 * @param[in] info   a line in a GTF file with information about feature (gene,
 *                   exon, transcript, etc.), start, and end in appropriate
 *                   locations. Ninth field "attribute" should contain
 *                   "transcript_id: " followed by transcript ID for optional
 *                   transcriptome-order output.
 *
 * @return New Sequence struct populated with available information. Returns
 * NULL if info does not contain the expected number of fields or if numerical
 * fields start and end cannot be parsed.
 */
Sequence get_sequence(string info) {
    Sequence seq;
    seq.start = -1;
    
    vector<string> fields = parse_tsv(info);
    if (fields.size() < NUM_GTF_ELT) {
        return seq;
    }
    
    seq.seqname = fields[0];
#if USING_GTF
    seq.feature = fields[2];
#else
    seq.feature = "exon";
#endif
    
#if USING_GTF
    seq.id = "";
    string id_start = ID_START;
    string id_end = ID_END;
    if (id_end.size() == 1) {
        int start = fields[8].find(id_start) + id_start.length();
        int end = fields[8].find(id_end, start + 1);
        if (start != string::npos && end != string::npos) {
            seq.id = fields[8].substr(start, end - start);
        }
    }
#else
    seq.id = fields[8];
#endif
    
    try {
        /* Using 0-indexed, half-interval in compliance with seqan. */
        seq.start = stoi(fields[3]) - 1;
        seq.end = stoi(fields[4]);
    }
    catch (invalid_argument &e) {
        seq.start = -1;
    }
    return seq;
}

/**
 * @brief Reads annotated sequences in (GTF) file. Attribute (9th) field should
 * contain "transcript_id: " followed by the transcript ID with no spaces.
 *
 * Note: assumes that information about a single chromosome/scaffold will not
 * be spread across multiple GTF files.
 *
 * @param file                  name of (GTF) file containing annotated
 *                              sequences
 * @param transcripts           (empty) vector to fill with sequence information
 * @param verbose               if 1, output warning messages to std::cerr
 * @param seq_count_start       number at which to start indexing transcripts.
 *                              if this is the first GTF, it should be 0, else
 *                              wherever the last GTF stopped
 *
 * @return                      1 if file fails to open, else 0
 */
int readGFF(string file, unordered_map<uint64_t, uint64_t> index_map,
            vector<vector<Exon>*> &exons, uint64_t &transcript_count,
            int verbose) {
     
    cout << "  Reading " << file << "..." << flush;
    
    string prev_ref = "", prev_transcript_id = "";
    unordered_map<string, Exon> chrom;
    
    uint64_t line_count = 0;
    uint64_t lines = get_line_count(file);
    if (lines == -1) {
        return 1;
    }
    uint64_t ten_percent = lines / 10;
    
    seqan::GffFileIn gff;
    if (!seqan::open(gff, file.c_str())) {
        return 1;
    }
    seqan::GffRecord rec;
    while(!seqan::atEnd(gff)) {
        ++line_count;
        seqan::readRecord(rec, gff);
#if 0
        if (verbose && lines > MIN_UPDATE && line_count % ten_percent == 0) {
            cout << "    " << line_count / ten_percent << "0% done" << endl;
        }
#endif   
        if (string(seqan::toCString(rec.ref)).compare(".")  == 0) {
            if (verbose) {
                cerr << endl << "    WARNING: Line " << line_count;
                cerr << " contains no information in seqid field. It will not ";
                cerr << "be used.";
            }
        }
        else if (lower(seqan::toCString(rec.type)).compare("exon") == 0) {
            string key = to_string(rec.beginPos) + "," + to_string(rec.endPos);
            string transcript_id = "";
            for (int i = 0; i < length(rec.tagNames); ++i) {
                if (lower(seqan::toCString(rec.tagNames[i])).compare(
                            "transcript_id") == 0) {
                    transcript_id = lower(seqan::toCString(rec.tagValues[i]));
                }
            }
            if (transcript_id.size() == 0) {
                cerr << endl << "    WARNING: Could not find transcript_id tag";
                cerr << " for line " << line_count << ".";
            }

            if (prev_ref.compare(seqan::toCString(rec.ref)) == 0) {
                if (prev_transcript_id.compare(transcript_id) != 0) {
                    ++transcript_count;
                    prev_transcript_id = transcript_id;
                }
                
                uint64_t id;
                if (!index_map.empty()) {
                    try {
                        id = index_map.at(transcript_count);
                    }
                    catch (out_of_range &err) {
                        id = -1;
                        cerr << endl << "    WARNING: unable to map GFF ";
                        cerr << "transcript " << transcript_count;
                        cerr << " to transcriptome. It will show up in EC as ";
                        cerr << id <<  ".";
                    }
                }
                else {
                    id = transcript_count;
                }
                
                if (chrom.find(key) == chrom.end()) {
                    chrom.emplace(key, Exon(rec));
                }
                chrom.at(key).transcripts->push_back(id);
            } else {
                prev_ref = seqan::toCString(rec.ref);
                if (!chrom.empty()) {
                    exons.push_back(map_values(chrom));
                    chrom.clear();
                }
                ++transcript_count;
                prev_transcript_id = transcript_id;
                
                uint64_t id;
                if (!index_map.empty()) {
                    try {
                        id = index_map.at(transcript_count);
                    }
                    catch (out_of_range &err) {
                        id = -1;
                        cerr << endl << "    WARNING: unable to map GFF ";
                        cerr << "transcript " << transcript_count;
                        cerr << " to transcriptome. It will show up in EC as ";
                        cerr << id <<  ".";
                    }
                }
                else {
                    id = transcript_count;
                }
                chrom.emplace(key, Exon(rec));
                chrom.at(key).transcripts->push_back(id);
            }
        }
    }
    
    // Add the exons from the last-encountered chromosome/scaffold.
    if (!chrom.empty()) {
        exons.push_back(map_values(chrom));
    }
    
    return 0;
}

int readGFFs(vector<string> &files, vector<string> &transcriptome,
             vector<vector<Exon>*> &exons, int verbose) {
    // 0-index the transcript counts. This will make count start at 0.
    uint64_t transcript_count = -1;
    
    // Map from this program's transcript indexing to kallisto's.
    unordered_map<uint64_t, uint64_t> *index_map
    = new unordered_map<uint64_t, uint64_t>;
    
    if (transcriptome.size() != 0) {
        int ret = get_index_to_kallisto_index(files, transcriptome, *index_map,
                                              verbose);
        if (ret == 1) {
            // It prints its own error message, so just return 1.
            return 1;
        }
    }
    
    for (uint i = 0; i < files.size(); ++i) {
        int ret = readGFF(files[i], *index_map, exons,
                          transcript_count, verbose);
        if (ret == 1) {
            cerr << "  ERROR: could not read " << files[i] << endl;
            return 1;
        }
    }
    
    delete index_map;
    
    return transcript_count;
}
