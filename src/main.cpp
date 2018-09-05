#include <iostream>
#include <fstream>
#include <time.h>
#include <getopt.h>
#include <seqan/bam_io.h>
#include <seqan/gff_io.h>
#include "TCC_Matrix.hpp"
#include "Mapper.hpp"
#include "common.hpp"
using namespace std;

/**
 * Given a GFF file, checks that 1) it will not needlessly splice transcripts
 * into exons, 2) it is grouped by chromosome, 3) it is grouped by
 * transcript id, and 4) all exons are in order (within transcript, by strand).
 *
 * Also outputs the max number of exons per transcript.
 *
 * NOTE: only (2) and (3) do not rely on others to be true. All others rely on
 * (3).
 *
 * @param gtf       Query GFF file.
 * @param flag      0 to check all, n to stop as soon as requirement n failed.
 * @return          -1 if error occurs, 0 if GFF file is proper, 1 if not.
 */
int checkGFF(string gtf, int flag, bool verbose, string outfile) {
    ifstream in(gtf);
    if (!in.is_open()) { return -1; }
    ofstream out;
    if (outfile.size() != 0) {
        out.open(outfile);
        if (!out.is_open()) { return -1; }
    }

    string input;
    int line = 0;
    unordered_set<string> chromosomes, transcripts;
    string curr_chrom, curr_transcript_id;
    int prev_exon_end, exons = 0, exon_entries = 0;
    bool strand; /* True if +, false if -. */
    bool properly_spliced = true, grouped_chrom = true,
         grouped_transcript = true, exons_in_order = true;
    int max_exons_per_transcript = 0;
    vector<int> exon_counts;

    while (getline(in, input)) {
        ++line;
        if (input.size() == 0 || input[0] == '#') {
            continue;
        }

        vector<string> inp = parseString(input, "\t", 0);
        
        /* inp[2] is the feature type, i.e. exon, gene, transcript. */
        if (inp[2].compare("exon") != 0) {
            continue;
        }

        ++exon_entries;

        /* inp[0] is the name of the chromosome/scaffold. */
        if (curr_chrom.compare(inp[0]) != 0) {
            if (chromosomes.find(inp[0]) == chromosomes.end()) {
                chromosomes.emplace(inp[0]);
                curr_chrom = inp[0];
            } else {
                grouped_chrom = false;
                if (flag == 2) {
                    break;
                }
            }
        }

        /* inp[8] is the `attributes` field. */
        string transcript_id;
        string tofind = "transcript_id \"";
        int start = inp[8].find(tofind);
        if (start == string::npos) {
            cerr << "Unable to find transcript ID on line " << line << endl;
            return -1;
        }
        start += tofind.size();
        int end = inp[8].find("\";", start);
        transcript_id = inp[8].substr(start, end - start);

        if (curr_transcript_id.compare(transcript_id) != 0) {
            if (transcripts.find(transcript_id) == transcripts.end()) {
                transcripts.emplace(transcript_id);
                curr_transcript_id = transcript_id;
                max_exons_per_transcript = max(max_exons_per_transcript, exons);
                exon_counts.push_back(exons);
                exons = 1;
                /* inp[6] is the strand. */
                if (inp[6].size() != 0 && inp[6][0] == '+') {
                    strand = true;
                } else if (inp[6].size() != 0 && inp[6][0] == '-') {
                    strand = false;
                } else {
                    cerr << "Unable to parse " << inp[6] << " on line " << line;
                    cerr << ". Expected either + or -." << endl;
                    return -1;
                }
                /* inp[3] is the start of the entry, and inp[4] the end. */
                if (strand && isNumber(inp[4])) {
                    prev_exon_end = stoi(inp[4]);
                } else if (!strand && isNumber(inp[3])) {
                    prev_exon_end = stoi(inp[3]);
                } else {
                    cerr << "Unable to parse " << inp[3] << " or " << inp[4];
                    cerr  << " on line " << line << ". Expected an int.";
                    cerr << endl;
                    return -1;
                }
            } else {
                grouped_transcript = false;
                if (flag == 3) {
                    break;
                }
            }
        }

        /* inp[3] is the start of the entry, and inp[4] the end. */
        else {
            ++exons;

            int exon_start;
            if (strand && isNumber(inp[3])) {
                exon_start = stoi(inp[3]);
            } else if (!strand && isNumber(inp[4])) {
                exon_start = stoi(inp[4]);
            } else {
                cerr << "Unable to parse " << inp[3] << " or " << inp[4];
                cerr  << " on line " << line << ". Expected an int.";
                cerr << endl;
                return -1;
            }

            if (prev_exon_end == exon_start) {
                properly_spliced = false;
                if (flag == 1) {
                    break;
                }
            } else if ((strand && prev_exon_end > exon_start) 
                    || (!strand && prev_exon_end < exon_start)) {
                exons_in_order = false;
                if (flag == 4) {
                    break;
                }
            } else if (strand && isNumber(inp[4])) {
                prev_exon_end = stoi(inp[4]);
            } else if (!strand && isNumber(inp[3])) {
                prev_exon_end = stoi(inp[3]);
            } else {
                cerr << "Unable to parse " << inp[3] << " or " << inp[4];
                cerr  << " on line " << line << ". Expected an int.";
                cerr << endl;
                return -1;
            }
        }
    }

    if (verbose ) {
        if (getline(in, input)) {
            cout << "Stopped at line " << line << " (before end of file)." << endl;
        } else {
            cout << "At end of file." << endl;
        }

        cout << "No contiguous exons: " << flush;
        if (properly_spliced) cout << "PASSED" << endl;
        else cout << "FAILED" << endl;
        cout << "Grouped by chromosome: " << flush;
        if (grouped_chrom) cout << "PASSED" << endl;
        else cout << "FAILED" << endl;
        cout << "Grouped by transcript: " << flush;
        if (grouped_transcript) cout << "PASSED" << endl;
        else cout << "FAILED" << endl;
        cout << "Exons in order: " << flush;
        if (exons_in_order) cout << "PASSED" << endl;
        else cout << "FAILED" << endl;
        cout << "Max number of exons per transcript: " << max_exons_per_transcript;
        cout << endl;
        cout << "Total number of exon entries: " << exon_entries << endl;
    }

    if (outfile.size() != 0) {
        for (int i = 1; i < exon_counts.size(); ++i) {
            out << exon_counts[i] << endl;
        }
    }

    return properly_spliced && grouped_chrom && grouped_transcript
        && exons_in_order;
}

/**
 * @brief Tests whether a file can be opened. Can test for read only, write
 * only, and read/write files.
 *
 * @param filename          name of file to open
 * @param mode              0: read only, 1: write only, 2: read/write
 * @return                  true if file can be opened, else false. if input
 *                          of parameter mode was invalid, returns false.
 */
bool testOpen(string filename, int mode) {
    if (mode == 0) {
        ifstream f(filename);
        if (f.is_open()) {
            f.close();
            return true;
        }
        else { return false; }
    }
    else if (mode == 1) {
        ofstream f(filename, ofstream::app);
        if (f.is_open()) {
            f.close();
            return true;
        }
        else { return false; }
    }
    else if (mode == 2) {
        fstream f(filename, ofstream::app);
        if (f.is_open()) {
            f.close();
            return true;
        }
        else { return false; }
    }
    else { return false; }
}

void printTime(time_t time) {
    string s = to_string((int)(time / 3600));
    if (s.size() == 1) { cout << "0"; }
    cout << s << ":";
    s = to_string((int)((time % 3600) / 60));
    if (s.size() == 1) { cout << "0"; }
    cout << s << ":";
    s = to_string((int)((time % 3600) % 60));
    if (s.size() == 1) { cout << "0"; }
    cout << s << endl;
}

void usage() {
    cerr << "Usage: thing [options]* -g <GFF> -S <BAM/SAM> [-o output]" << endl
    << "  <GFF>                     Comma-separated list of GFFs." << endl
    << "  <SAM/BAM>                 Comma-separated list of SAM/BAM files "
    << "sorted by genomic coordinate." << endl
    << "  <output>                  Prefix of output files (defaults to "
    << "`matrix`)" << endl
    << endl << "Options:" << endl
    << "  -U                        Indicate that reads are unpaired." << endl
    << "  -t, --transcriptome <fa>  Change transcript indexing to match that "
    << "of these comma-separated transcriptome files. Use to directly compare "
    << "output to that of kalisto pseudo." << endl
    << "  -e <EC>                   Use the ECs in this file. Any additional "
    << "ECs will be appended." << endl
    << "  -p <threads>              Number of threads to use. Currently not "
    << "operational." << endl
    << "  --full-matrix             Output full (not sparse) matrix." << endl
    << "  -u, --unmapped <SAM/BAM>  Output unmapped reads to files <SAM/BAM>."
    << " Must provide one for each input SAM/BAM file."
    << endl
    << "  --check-gff               Check GFFs only." << endl
    << endl;
}

int main(int argc, char **argv) {
    time_t startTime = time(0);
    vector<string> gff, bam, fa, unmapped;
    string outprefix = "matrix", ec = "";
    bool paired = true, full = false, checkGFFOnly = false,
         pgProvided = false, genomebam = false, rapmap = false;
    int threads = 1;
    
    /* Parse options. */
    struct option opts[] = {
        {"GFF", required_argument, 0, 'g'},
        {"SAM/BAM", required_argument, 0, 'S'},
        {"output", required_argument, 0, 'o'},
        {"Unpaired", no_argument, no_argument, 'U'},
        {"transcriptome", required_argument, 0, 't'},
        {"EC", required_argument, 0, 'e'},
        {"threads", required_argument, 0, 'p'},
        {"full-matrix", no_argument, no_argument, 'f'},
        {"unmapped", required_argument, 0, 'u'},
        {"check-gff", no_argument, no_argument, 'G'},
        {"kallisto", no_argument, no_argument, 'k'},
        {"RapMap", no_argument, no_argument, 'R'}
    };
    int opt_index = 0;
    while (true) {
        int c = getopt_long(argc, argv, "g:S:o:Ut:p:e:fu:kR", opts, &opt_index);
        if (c == -1) { break; }
        switch (c) {
            case 'g':   gff = parseString(optarg, ",", 0); break;
            case 'S':   bam = parseString(optarg, ",", 0); break;
            case 'o':   outprefix = optarg; break;
            case 'U':   paired = false; break;
            case 't':   fa = parseString(optarg, ",", 0); break;
            case 'p':   threads = atoi(optarg); break;
            case 'e':   ec = optarg; break;
            case 'f':   full = true; break;
            case 'u':   unmapped = parseString(optarg, ",", 0); break;
            case 'G':   checkGFFOnly = true; break;
            case 'k':   genomebam = true; pgProvided = true; break;
            case 'R':   rapmap = true; pgProvided = true; break;
        }
    }

    /* Make sure all required files are present. */
    if ((checkGFFOnly && gff.size() == 0)
           || (!checkGFFOnly && bam.size() == 0)
           || (unmapped.size() != 0 && unmapped.size() != bam.size())) {
        usage();
        return 1;
    }

    /* Test-open files. */
    for (auto file = gff.begin(); file != gff.end(); ++file) {
        seqan::GffFileIn f;
        if (!seqan::open(f, file->c_str())) {
            cerr << "ERROR: failed to open GFF file " << *file << endl;
            return 1;
        }
        seqan::GffRecord rec;
        seqan::readRecord(rec, f);
        seqan::close(f);
    }
    for (auto file = bam.begin(); file != bam.end(); ++file) {
        seqan::BamFileIn f;
        if (!seqan::open(f, file->c_str())) {
            cerr << "ERROR: failed to open SAM/BAM file " << *file << endl;
            return 1;
        }
        seqan::BamHeader head;
        seqan::readHeader(head, f);
        seqan::BamAlignmentRecord rec;
        seqan::readRecord(rec, f);
        seqan::close(f);
    }
    for (auto file = fa.begin(); file != fa.end(); ++file) {
        if (!testOpen(*file, 0)) {
            cerr << "ERROR: failed to open FASTA file " << *file << endl;
            return 1;
        }
    }
    if (ec.size() != 0 && !testOpen(ec, 0)) {
            cerr << "ERROR: failed to open reference EC " << ec << endl;
            return 1;
    }
    if (!testOpen(outprefix + ".ec", 1)) {
        cerr << "ERROR: failed to open output file " << outprefix << ".ec"
            << endl;
        return 1;
    }
    if (!testOpen(outprefix + ".tsv", 1)) {
        cerr << "ERROR: failed to open output file " << outprefix << ".tsv"
            << endl;
        return 1;
    }
    if (!testOpen(outprefix + ".cells", 1)) {
        cerr << "ERROR: failed to open output file " << outprefix << ".cells"
            << endl;
        return 1;
    }
    for (auto file = unmapped.begin(); file != unmapped.end(); ++file) {
        if (!testOpen(*file, 1)) {
            cerr << "ERROR: failed to open output file " << *file << endl;
            return 1;
        }
    }

    /* Check GFF formatting. */
    cout << "Checking GFF formatting..." << endl;
    for (auto file = gff.begin(); file != gff.end(); ++file) {
        if (!checkGFF(*file, 0, checkGFFOnly, ""))
        {
            cerr << *file << " incorrectly formatted." << endl;
            return 1;
        }
    }
    if (checkGFFOnly) { return 0; }

    /* Map and write */
    Mapper mapper(gff, bam, fa, paired, unmapped.size() != 0);
    cout << "Mapping reads..." << endl;
    mapper.mapReads(threads, pgProvided, genomebam, rapmap);
    cout << "Writing to file..." << endl;
    mapper.writeToFile(outprefix, unmapped, full, ec);
    
    printTime(time(0) - startTime);
    return 0;
}

