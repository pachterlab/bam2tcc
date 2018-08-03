/**
 * Functions for reading GFF files.
 */
#include <seqan/gff_io.h>
#include "gff_io.hpp"
#include "util.hpp"
#include "kallisto_util.hpp"
using namespace std;

/**
 * @brief Returns a vector containing all exons (values) in map.
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
 * @brief Reads a GFF file and fills `exons` appropriately.
 *
 * Note: assumes that information about a single chromosome/scaffold will not
 * be spread across multiple GFF files.
 *
 * @param file                  Name of GFF file containing annotated sequences.
 *
 * @param index_map             A map from transcripts' indices in GFF file to
 * transcripts' indices in transcriptome file. See get_index_to_kallisto_index
 * in kallisto_util.cpp. Optional.
 *
 * @param exons                 Vector to be filled with information in GFF.
 *
 * @param seq_count_start       Number at which to start indexing transcripts.
 * If this is the first GFF, it should be 0, else wherever the last GFF stopped.
 *
 * @param verbose               If 1, output warning messages to std::cerr.
 *
 * @return                      -1 if file fails to open, else 0.
 */
int readGFF(string file, unordered_map<int, int> index_map,
            unordered_map<string, vector<Exon>*> &exons, int &transcript_count,
            int verbose) {
     
    cout << "  Reading " << file << "..." << flush;
    
    string prev_ref = "", prev_transcript_id = "";
    unordered_map<string, Exon> chrom;
    
    int line_count = 0;
    
    seqan::GffFileIn gff;
    if (!seqan::open(gff, file.c_str())) {
        return -1;
    }
    seqan::GffRecord rec;
    while(!seqan::atEnd(gff)) {
        ++line_count;
        seqan::readRecord(rec, gff);
        
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
                
                int id;
                if (!index_map.empty()) {
                    try {
                        id = index_map.at(transcript_count);
                    }
                    catch (out_of_range &err) {
                        /* This shouldn't happen. If it does, it indicates a bug
                         * in the get_index_to_kallisto_index function in
                         * kallisto_util.cpp. */
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
                if (!chrom.empty()) {
                    exons.emplace(lower(prev_ref), map_values(chrom));
                    chrom.clear();
                }
                prev_ref = seqan::toCString(rec.ref);
                ++transcript_count;
                prev_transcript_id = transcript_id;
                
                int id;
                if (!index_map.empty()) {
                    try {
                        id = index_map.at(transcript_count);
                    }
                    catch (out_of_range &err) {
                        /* This shouldn't happen. If it does, it indicates a bug
                         * in the get_index_to_kallisto_index function in
                         * kallisto_util.cpp. */
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
    
    /* Add the exons from the last-encountered chromosome/scaffold. */
    if (!chrom.empty()) {
        exons.emplace(lower(prev_ref), map_values(chrom));
    }
 
    return 0;
}

/**
 * @brief Read in all GFF files named in vector `files`.
 *
 * @param files         Vector containing the names of query GFF files.
 *
 * @param transcriptome Vector containing the names of query transcriptome
 * files. If non-empty, readGFFs will use the transcript indexing in the
 * transcriptome files. If multiple files are provided, the first transcript
 * of the first file is index 0, and the indexing of remainng files start
 * where the previous one left off.
 *
 * @param exons                 Vector to be filled with information in GFFs.
 *
 * @param verbose               If 1, output warning messages to std::cerr.
 *
 * @return                      -1 if error occurs, else number of transcripts
 *                              read.
 */
int readGFFs(vector<string> &files, vector<string> &transcriptome,
             unordered_map<string, vector<Exon>*> &exons, int verbose) {

    cout << "Reading GFFs..." << endl;

    /* 0-index the transcript counts. This will make count start at 0. */
    int transcript_count = -1;
    
    /* Map from this program's transcript indexing to kallisto's. */
    unordered_map<int, int> *index_map = new unordered_map<int, int>;
    
    if (transcriptome.size() != 0) {
        int ret = get_index_to_kallisto_index(files, transcriptome, *index_map,
                                              verbose);
        if (ret == -1) {
            /* It prints its own error message, so just return 1. */
            return -1;
        }
    }
    
    for (int i = 0; i < files.size(); ++i) {
        int ret = readGFF(files[i], *index_map, exons,
                          transcript_count, verbose);
        if (ret == -1) {
            cerr << "  ERROR: could not read " << files[i] << endl;
            return -1;
        }
    }
    
    delete index_map;
    ++transcript_count;
    cout << "  read " << transcript_count << " transcripts" << endl;
    return transcript_count;
}

