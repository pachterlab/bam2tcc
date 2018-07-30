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
int readGTF(string file, unordered_map<uint64_t, uint64_t> index_map,
            vector<vector<Exon>*> &exons, uint64_t &transcript_count,
            int verbose) {
    
    // Holds one line of the file. For use with getline.
    string inp;
    
    // seqname and id of previous exon. For use with seq_count.
    string prev_seqname = "", prev_id = "";
    
    // Maps "start,end" to exon in ONE chromosome or scaffold. Add info to
    // vector exons and clear whenever a new chromosome or scaffold is
    // encountered.
    unordered_map<string, Exon> chrom;
    
    // Number of current line being read. For use in error messages if a line
    // cannot be parsed.
    uint64_t line_count = 0;
    
    // Total number of lines in file. If verbose, prints update messages every
    // 10% done.
    uint64_t lines = get_line_count(file);
    if (lines == -1) { // File could not be opened! Return 1.
        return 1;
    }
    // 10% of the file. If verbose, prints update messages every ten_percent
    // lines.
    uint64_t ten_percent = lines / 10;
    
    /* Open GTF files. Return 1 if error occurs */
    ifstream f(file);
    if (!f.is_open()) {
        return 1;
    }
    
    cout << "  Reading " << file << "..." << endl;
    
    /* Start while loop to go through GTF file */
    while(getline(f, inp)) {
        /* Update line count and print update message */
        ++line_count;
        if (verbose && lines > MIN_UPDATE && line_count % ten_percent == 0) {
            cout << "    " << line_count / ten_percent << "0% done" << endl;
        }
        
        /* Skip header and blank lines */
        if (inp.size() == 0 || inp[0] == '#') {
            continue;
        }
        
        /* Parse input and place information in seq */
        inp = lower(inp);
        Sequence seq = get_sequence(inp);
        
        /* Look at values and update vector exons, map chrom, and various
         counters appropriately */
        if (seq.start == -1) {
            if (verbose) {
                cerr << "    WARNING: Failed to read line ";
                cerr << line_count << endl;
            }
        }
        else if (seq.seqname.size() == 0) {
            if (verbose) {
                cerr << "    WARNING: Line " << line_count;
                cerr << " contains no information in seqname field" << endl;
            }
        }
        else if (seq.feature.compare("exon") == 0) {
            // Key for map chrom.
            string key = to_string(seq.start) + "," + to_string(seq.end);
            
            // Are we looking at the same chromosome/scaffold? If so, we update
            // the relevant exon in map chrom.
            if (prev_seqname.compare(seq.seqname) == 0) {
                // Only update seq_count if we're looking at an exon in
                // a new transcript.
                if (prev_id.compare(seq.id) != 0) {
                    ++transcript_count;
                    prev_id = seq.id;
                }
                
                uint64_t id;
                if (!index_map.empty()) {
                    try {
                        id = index_map.at(transcript_count);
                    }
                    catch (out_of_range &err) {
                        id = -1;
                        cerr << "  WARNING: unable to map index ";
                        cerr << transcript_count << " to transcriptome's. ";
                        cerr << "It will show up in EC as " << id <<  ".";
                        cerr << endl;
                    }
                }
                else {
                    id = transcript_count;
                }
                
                if (chrom.find(key) == chrom.end()) {
                    chrom.emplace(key, Exon(seq));
                }
                chrom.at(key).transcripts->push_back(id);
            }
            
            // It's a new chromosome/scaffold, so we need to clear map chrom
            // and create a new vector in which to place the exon.
            else {
                // Update seqname.
                prev_seqname = seq.seqname;
                
                // If we have something in chrom, we add it to vector exons
                // and clear it.
                if (!chrom.empty()) {
                    exons.push_back(map_values(chrom));
                    chrom.clear();
                }
                
                // This is definitely a new transcript, so update.
                ++transcript_count;
                prev_id = seq.id;
                
                uint64_t id;
                if (!index_map.empty()) {
                    try {
                        id = index_map.at(transcript_count);
                    }
                    catch (out_of_range &err) {
                        id = -1;
                        cerr << "  WARNING: unable to map index ";
                        cerr << transcript_count << " to transcriptome's. ";
                        cerr << "It will show up in EC as " << id << ".";
                        cerr << endl;
                    }
                }
                else {
                    id = transcript_count;
                }
                
                // Add the new exon to the now-empty map chrom.
                chrom.emplace(key, Exon(seq));
                chrom.at(key).transcripts->push_back(id);
            }
        }
    }
    
    // Add the exons from the last-encountered chromosome/scaffold.
    if (!chrom.empty()) {
        exons.push_back(map_values(chrom));
    }
    
    f.close();
    return 0;
}

/**
 *
 */
int readGTFs(vector<string> &files, vector<string> &transcriptome,
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
        int ret = readGTF(files[i], *index_map, exons,
                          transcript_count, verbose);
        if (ret == 1) {
            cerr << "  ERROR: could not read " << files[i] << endl;
            return 1;
        }
    }
    
    delete index_map;
    
    return transcript_count;
}
