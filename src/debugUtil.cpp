#include <fstream>
#include <iostream>
#include <unordered_map>
#include <unordered_set>
#include <set>

#include "common.hpp"
#include "FileUtil.hpp"
#include "TCC_Matrix.hpp"
using namespace std;

#define TRANSCRIPT_NAME_START_CHAR '>'
#define TRANSCRIPT_NAME_END_CHAR ' '

int no_zeros(string inname, string outname) {
    ifstream inec(inname + ".ec");
    ifstream intsv(inname + ".tsv");
    if (!inec.is_open() || !intsv.is_open()) {
        cerr << "ERROR: cannot open input files" << endl;
        return 1;
    }
    ofstream outec(outname + ".ec");
    ofstream outtsv(outname + ".tsv");
    if (!outec.is_open() || !outtsv.is_open()) {
        cerr << "ERROR: cannot open output files" << endl;
        return 1;
    }
    
    string inpec, inptsv;
    while (getline(inec, inpec) && getline(intsv, inptsv)) {
        vector<string> tcc = parseString(inptsv, "\t", 0);
        int flag = 0;
        for (uint i = 1; i < tcc.size(); ++i) {
            if (tcc[i].compare("0") != 0) {
                flag = 1;
                break;
            }
        }
        if (flag) {
            outec << inpec << endl;
            outtsv << inptsv << endl;
        }
    }
    
    inec.close();
    intsv.close();
    outec.close();
    outtsv.close();
    
    return 0;
}

/**
 * This, uh, turns the kallisto single-end output full matrix into a sparse
 * matrix. Only works when there's only one column. lol
 */
int moar_zeroes(string infile, string outfile) {
    ifstream in(infile);
    if (!in.is_open()) {
        cerr << infile << "  failed to open" << endl;
        return 1;
    }
    ofstream out(outfile);
    if (!out.is_open()) {
        cerr << outfile << " failed to open" << endl;
        in.close();
        return 1;
    }

    string inp;
    while (getline(in, inp)) {
        vector<string> tsv = parseString(inp, "\t", 0);
        if (stoi(tsv[1]) != 0) {
            out << tsv[0] << "\t0\t" << tsv[1] << endl;
        }
    }

    in.close();
    out.close();
    return 0;
}

int count(string inname) {
    ifstream in(inname + ".tsv");
    if (!in.is_open()) {
        cerr << "ERROR: unable to open " << inname << ".tsv" << endl;
        return 1;
    }
    
    string inp;
    getline(in, inp);
    vector<string> tcc = parseString(inp, "\t", 0);
    
    vector<uint> count;
    for (uint i = 0; i < tcc.size(); ++i) {
        if (isNumber(tcc[i])) {
            count.push_back((uint) stoi(tcc[i]));
        } else {
            count.push_back(0);
        }
    }
    
    while (getline(in, inp)) {
        tcc = parseString(inp, "\t", 0);
        for (uint i = 0; i < tcc.size(); ++i) {
            if (isNumber(tcc[i])) {
                count[i] += stoi(tcc[i]);
            }
        }
    }
    
    for (uint i = 0; i < count.size(); ++i) {
        cout << count[i] << '\t';
    }
    cout << endl;
    
    in.close();
    
    return 0;
}

int intersection(string inname1, string inname2, string outname) {
    ifstream in1(inname1 + ".tsv");
    if (!in1.is_open()) {
        cerr << "ERROR: unable to open" << inname1 << ".tsv" << endl;
        return 1;
    }
    
    ifstream in2(inname2 + ".tsv");
    if (!in2.is_open()) {
        cerr << "ERROR: unable to open" << inname2 << ".tsv" << endl;
        return 1;
    }
    
    ofstream out(outname + ".tsv");
    if (!out.is_open()) {
        cerr << "ERROR: unable to open" << outname << ".tsv" << endl;
        return 1;
    }
    
    string inp1, inp2;
    if (!getline(in1, inp1)) {
        in1.close();
        in2.close();
        out.close();
        return 0;
    }
    if (!getline(in2, inp2)) {
        in1.close();
        in2.close();
        out.close();
        return 0;
    }
    int index1 = stoi(parseString(inp1, "\t", 1)[0]);
    int index2 = stoi(parseString(inp2, "\t", 1)[0]);
    while (true) {
        if (index1 == index2) {
            out << inp1 << endl << inp2 << endl << endl;
            if (!getline(in1, inp1)) {
                break;
            }
            index1 = stoi(parseString(inp1, "\t", 1)[0]);
            if (!getline(in2, inp2)) {
                break;
            }
            index2 = stoi(parseString(inp2, "\t", 1)[0]);
        }
        else if (index1 < index2) {
            if (!getline(in1, inp1)) {
                break;
            }
            index1 = stoi(parseString(inp1, "\t", 1)[0]);
        }
        else {
            if (!getline(in2, inp2)) {
                break;
            }
            index2 = stoi(parseString(inp2, "\t", 1)[0]);
        }
    }
    
    in1.close();
    in2.close();
    out.close();
    
    return 0;
}

int stats(string inname1, string inname2) {
    ifstream in1(inname1 + ".tsv");
    ifstream in2(inname2 + ".tsv");
    if (!in1.is_open() || !in2.is_open()) {
        cerr << "ERROR: cannot open input files" << endl;
        return 1;
    }
    
    int same = 0;
    int ten_per = 0;
    int fifty_per = 0;
    int wrong = 0;
    int false_zero1 = 0;
    int false_zero2 = 0;
    int total_nonzero = 0;
    int total = 0;
    string inp1, inp2;
    while (getline(in1, inp1) && getline(in2, inp2)) {
        int count1 = stoi(parseString(inp1, "\t", 2)[1]);
        int count2 = stoi(parseString(inp2, "\t", 2)[1]);
        ++total;
        if (count1 != 0 && count2 != 0) {
            ++total_nonzero;
            if (count1 == count2) {
                ++same;
            }
            else if (2. * abs(count1 - count2) / (count1 + count2) < .1) {
                ++ten_per;
            }
            else if (2. * abs(count1 - count2) / (count1 + count2) < .5) {
                ++fifty_per;
            }
            else if (2. * abs(count1 - count2) / (count1 + count2) >= 1) {
                ++wrong;
            }
        }
        if (count1 == 0 && count2 != 0) {
            ++false_zero1;
        }
        if (count1 != 0 && count2 == 0) {
            ++false_zero2;
        }
    }
    
    cout << "Out of " << total << " entries: " << endl;
    cout << "Same (non-zero): " << 100. * same / total_nonzero << "%" << endl;
    cout << "Within ten percent: " << 100. * ten_per / total_nonzero << "%" << endl;
    cout << "Within fifty percent: " << 100. * fifty_per / total_nonzero << "%" << endl;
    cout << "Wrong: " << 100. * wrong / total_nonzero << "%" << endl;
    cout << "False zeroes in " << inname1 << ": " << false_zero1 << endl;
    cout << "False zeroes in " << inname2 << ": " << false_zero2 << endl;
    
    in1.close();
    in2.close();
    return 0;
}

int transcript_ids(string infile, string outfile, string transcriptome) {
    ifstream in(infile);
    if (!in.is_open()) {
        cout << "File " << infile << " could not be opened." << endl;
        return -1;
    }
    ifstream trans(transcriptome);
    if (!trans.is_open()) {
        cout << "Transcriptome " << transcriptome << " could not be opened.";
        cout << endl;
        return -1;
    }
    ofstream out(outfile);
    if (!out.is_open()) {
        cout << "File " << outfile << " could not be opened." << endl;
        return -1;
    }
    
    unordered_map<string, string> *map = new unordered_map<string, string>;
    string inp;
    uint64_t count = 0;
    while (getline(trans, inp)) {
        if (inp.size() == 0 || inp[0] != '>') {
            continue;
        }
        
        inp = lower(inp);
        uint first_space = inp.find(' ');
        // Exclude the first character ('>') and the END CHAR.
        map->emplace(to_string(count), inp.substr(1, first_space - 1));
        ++count;
    }
    
    while (getline(in, inp)) {
        inp = lower(inp);
        vector<string> t1 = parseString(inp, "\t", 0);
        out << t1[0] << '\t';
        vector<string> t2 = parseString(t1[1], ",", 0);
        for (uint i = 0; i < t2.size(); ++i) {
            try {
                out << map->at(t2[i]) << ',';
            }
            catch (out_of_range &e) {
                out << "-1,";
            }
        }
        out << endl;
    }
    
    in.close();
    trans.close();
    out.close();
    delete map;
    return 0;
}

int find_in_transcriptome_transcript(string tofind, int transcript_num,
                                     string infile) {
    ifstream in(infile);
    if (!in.is_open()) {
        cout << "File " << infile << " could not be opened." << endl;
        return -1;
    }
    int count = -1;
    string inp;
    while (getline(in, inp)) {
        if (inp.size() != 0 && inp[0] == '>') {
            ++count;
        }
        if (count == transcript_num) {
            break;
        }
    }
    string seq = "";
    while (getline(in, inp)) {
        if (inp.size() != 0 && inp[0] == '>') {
            break;
        }
        seq += inp;
    }
    
    if (seq.find(tofind) != string::npos) {
            cout << "Sequence found" << endl;
    }
    else {
        cout << "Sequence not found" << endl;
    }
    return 0;
}

int add_transcript_names(string infile, string outfile, string transcriptome) {
    ifstream in(infile);
    ofstream out(outfile);
    ifstream trans(transcriptome);
    if (!in.is_open() || !out.is_open() || !trans.is_open()) {
        cerr << "Could not open a file" << endl;
        return -1;
    }
    string inp;
    getline(in, inp);
    out << inp << "\ttranscript_name" << endl;
    while (getline(trans, inp)) {
        if (inp.size() == 0 || inp[0] != '>') {
            continue;
        } else {       
            int start = inp.find(TRANSCRIPT_NAME_START_CHAR) + 1;
            int end = inp.find(TRANSCRIPT_NAME_END_CHAR, start);
            string transcriptEC = "\t";
            if (start == string::npos || end == string::npos) {
                cerr << "  ERROR: unexpected input" << endl;
                transcriptEC += "N/A";
            } else {
                transcriptEC += inp.substr(start, end - start);
            }
            getline(in, inp);
            out << inp << transcriptEC << endl;
        }
    }
    while (getline(in, inp)) {
        out << inp << "\tN/A" << endl;
    }
    in.close();
    out.close();
    trans.close();
    return 0;
}

int fill_fastq(string sam, string fastq_in, string fastq_out) {
    ifstream in(sam);
    if (!in.is_open()) {
        cerr << "SAM file could not be opened" << endl;
        return 1;
    }
    ofstream out(fastq_out);
    if (!out.is_open()) {
        cerr << "fastq outfile could not be opened" << endl;
        return 1;
    }
    
    string inp;
    
    vector<string> reads;
    while (getline(in, inp)) {
        string r = inp.substr(0, inp.find('\t'));
        if (reads.size() == 0 || r.compare(reads[reads.size() - 1]) != 0) {
            reads.push_back(r);
        }
    }
    
    ifstream fastq(fastq_in);
    if (!fastq.is_open()) {
        cerr << "transcriptome could not be opened" << endl;
        return 1;
    }
    int found = 0;
    while (getline(fastq, inp)) {
        for (uint i = 0; i < reads.size(); ++i) {
            if (reads[i].compare(inp.substr(1, inp.find(' ') - 1)) == 0) {
                out << inp << endl;
                getline(fastq, inp);
                out << inp << endl;
                ++found;
                break;
            }
        }
        if (found == reads.size()) {
            break;
        }
    }
    if (found != reads.size()) {
        cout << reads.size() - found << " not found" << endl;
    }
    
    in.close();
    out.close();
    return 0;
}

int salmon_EQ_to_TCC(string infile, string outprefix, string transcriptome,
        string ref_ec) {
    ifstream in(infile);
    if (!in.is_open()) {
        cerr << "Unable to open " << infile << endl;
        return -1;
    }
    ifstream trans(transcriptome);
    if (!trans.is_open()) {
        cerr << "Unable to open " << transcriptome << endl;
        return -1;
    }
    ofstream cells(outprefix + ".cells");
    if (!cells.is_open()) {
        cerr << "Unable to open outfiles with prefix " << outprefix << endl;
        return -1;
    }
    cells << outprefix << endl;
    cells.close();

    string inp;
    int line_count = 0;
    unordered_map<string, int> *m1 = new unordered_map<string, int>;
    while(getline(trans, inp)) {
        if (inp.size() == 0 || inp[0] != '>') {
            continue;
        }
        inp = lower(inp);
        int end = inp.find(' ');
        if (end == string::npos) {
            cerr << "ERROR: unexpected input in transcriptome" << endl;
            return -1;
        }
        m1->emplace(inp.substr(1, end - 1), line_count);
        ++line_count;
    }

    getline(in, inp);
    int num_transcripts = stoi(inp);
    getline(in, inp);
    int num_ECs = stoi(inp);
    line_count = 0;
    unordered_map<int, int> *m = new unordered_map<int, int>;
    while (line_count < num_transcripts) {
        getline(in, inp);
        m->emplace(line_count, m1->at(lower(inp)));
        ++line_count;
    }

    line_count = 0;
    TCC_Matrix *matrix = new TCC_Matrix(1);
    while (line_count < num_ECs) {
        getline(in, inp);
        vector<string> EC = parseString(inp, "\t", 0);
        string stringEC = to_string(m->at(stoi(EC[1])));
        for (int i = 2; i < stoi(EC[0]) + 1; ++i) {
            stringEC += "," + to_string(m->at(stoi(EC[i])));
        }
        for (int i = 0; i < stoi(EC[EC.size() - 1]); ++i) {
            matrix->inc_TCC(stringEC, 0);
        }
        ++line_count;
    }
    in.close();

    vector<string> *kallisto_order = new vector<string>;
    unordered_set<string> *kallisto_ecs = new unordered_set<string>;
    int err = getECOrder(ref_ec, *kallisto_order,
                                    *kallisto_ecs);
    if (err != -1) {
        err = matrix->write_to_file_in_order_sparse(outprefix,
            *kallisto_order, *kallisto_ecs);
    }
    delete m;
    delete m1;
    delete matrix;
    delete kallisto_order;
    delete kallisto_ecs;
    return err;
}

int normal_transcriptome(string infile, string outfile) {
    ifstream in(infile);
    if (!in.is_open()) {
        cerr << "Unable to open " << infile << endl;
        return -1;
    }
    ofstream out(outfile);
    if (!out.is_open()) {
        cerr << "Unable to open " << outfile << endl;
        return -1;
    }
    string inp;
    while (getline(in, inp)) {
        if (inp.size() == 0 || inp[0] != '>') {
            out << inp << endl;
            continue;
        }
        int start = inp.find(' ');
        if (start == string::npos) {
            cerr << "Unexpected input in " << infile << ':' << endl;
            cerr << inp << endl;
            return -1;
        }
        out << '>' << inp.substr(start + 1, inp.size() - start + 1) << endl;
    }
    in.close();
    out.close();
    return 0;
}

int gene_table(string infile, string outfile) {
    ifstream in(infile);
    if (!in.is_open()) {
        cout << "Unable to open " << infile << endl;
        return -1;
    }
    ofstream out(outfile);
    if (!out.is_open()) {
        cout << "Unable to open " << outfile << endl;
        return -1;
    }

    string inp;
    while (getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '#') {
            continue;
        }

        vector<string> input = parseString(inp, "\t", 0);
        if (input.size() != 9) {
            cout << "Not enough fields in GFF" << endl;
            return -1;
        }
        if (input[2].compare("transcript") == 0) {
            string attribute = input[8];
            string find = "gene_id \"";
            int start = attribute.find(find);
            if (start == string::npos) {
                cout << "Gene ID not found:" << endl;
                cout << inp << endl;
                continue;
            } else {
                start += find.size();
            }
            int end = attribute.find("\";", start);
            out << attribute.substr(start, end - start) << '\t';
            find = "transcript_id \"";
            start = attribute.find(find);
            if (start == string::npos) {
                cout << "Transcript ID not found:" << endl;
                cout << inp << endl;
                continue;
            } else {
                start += find.size();
            }
            end = attribute.find("\";", start);
            out << attribute.substr(start, end - start) << endl;
        }
    }
    in.close();
    out.close();
    return 0;
}

int transcript_index(string transcriptome, string outfile) {
    ifstream in(transcriptome);
    if (!in.is_open()) {
        cout << "Unable to open " << transcriptome << endl;
        return -1;
    }
    ofstream out(outfile);
    if (!out.is_open()) {
        cout << "Unable to open " << outfile << endl;
        return -1;
    }
    int count = 0;
    string inp;
    while (getline(in, inp)) {
        if (inp.size() != 0 && inp[0] == '>') {
            int start = inp.find(TRANSCRIPT_NAME_START_CHAR);
            int end = inp.find(TRANSCRIPT_NAME_END_CHAR, start);
            if (start == string::npos || end == string::npos) {
                cerr << "  ERROR: unexpected input" << endl;
                out << count << '\t' << "N/A" << endl;
                ++count;
            }
            ++start;
            out << count << '\t' << inp.substr(start, end - start) << endl;
            ++count;
        }
    }
    in.close();
    out.close();
    return 0;
}

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
int checkGFF(string gtf, int flag=0, string outfile="") {
    ifstream in(gtf);
    if (!in.is_open()) {
        cerr << "Unable to open " << gtf << endl;
        return -1;
    }
    ofstream out;
    if (outfile.size() != 0) {
        out.open(outfile);
        if (!out.is_open()) {
            cerr << "Unable to open " << outfile << endl;
            return -1;
        }
    }

    string input;
    int line = 0;
    unordered_set<string> chromosomes;
    unordered_set<string> transcripts;
    string curr_chrom;
    string curr_transcript_id;
    int prev_exon_end;
    int exons = 0;
    int exon_entries = 0;
    bool strand; /* True if +, false if -. */
    bool properly_spliced = true;
    bool grouped_chrom = true;
    bool grouped_transcript = true;
    bool exons_in_order = true;
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

    if (outfile.size() != 0) {
        for (int i = 1; i < exon_counts.size(); ++i) {
            out << exon_counts[i] << endl;
        }
    }

    return properly_spliced && grouped_chrom && grouped_transcript
        && exons_in_order;
}

int timecourse1(string sam, string tsv, string outfile) {
    ifstream in(sam);
    if (!in.is_open()) {
        cout << "Failed to open " << sam << endl;
        return 1;
    }
    ifstream in2(tsv);
    if (!in2.is_open()) {
        cout << "Failed to open " << tsv << endl;
        return 1;
    }
    ofstream out(outfile, ofstream::app);
    if (!out.is_open()) {
        cout << "Failed to open " << outfile << endl;
        return 1;
    }

    string inp;
    string prevQName;
    int reads = 0;
    int alignments = 0;
    int alignedAlignments = 0;
    unordered_map<string, int> chromosomeCounts;
    while (getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '@') {
            continue;
        }
        vector<string> alignment = parseString(inp, "\t", 3);
        ++alignments;
        if (prevQName.compare(alignment[0])) {
            prevQName = alignment[0];
            ++reads;
        }
        if ((((uint) stoi(alignment[1])) & 0x4) == 0) {
            ++alignedAlignments;
        }
        if (chromosomeCounts.find(alignment[2]) == chromosomeCounts.end()) {
            chromosomeCounts.emplace(alignment[2], 1);
        } else {
            ++chromosomeCounts[alignment[2]];
        }
    }

    int mapped = 0;
    while (getline(in2, inp)) {
        mapped += stoi(parseString(inp, "\t", 3)[2]);
    }

    out << alignedAlignments << '\t' << alignments << '\t' << reads << '\t'
        << mapped << "\t" << chromosomeCounts.size() << '\t'
        << "[time]" << endl;

    in.close();
    in2.close();
    out.close();
    return 0;
}

int pull_times(string log) {
    ifstream in(log);
    if (!in.is_open()) {
        cout << "Unable to open " << log << endl;
        return 1;
    }
    string inp;
    while(getline(in, inp)) {
        vector<string> vec = parseString(inp, ":", 0);
        if (vec.size() != 3) {
            continue;
        }
        cout << 3600 * stoi(vec[0]) + 60 * stoi(vec[1]) + stoi(vec[2]) << endl;
    }
    in.close();
    return 0;
}

int timecourse(string sam, string tsv, string log, string outfile) {
    ifstream in(sam);
    if (!in.is_open()) {
        cout << "Failed to open " << sam << endl;
        return 1;
    }
    ifstream in2(tsv);
    if (!in2.is_open()) {
        cout << "Failed to open " << tsv << endl;
        return 1;
    }
    ifstream in3(log);
    if (!in3.is_open()) {
        cout << "Failed to open " << log << endl;
        return 1;
    }
    ofstream out(outfile, ofstream::app);
    if (!out.is_open()) {
        cout << "Failed to open " << outfile << endl;
        return 1;
    }

    string inp;
    string prevQName;
    int reads = 0;
    int alignments = 0;
    int alignedAlignments = 0;
    unordered_map<string, int> chromosomeCounts;
    while (getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '@') {
            continue;
        }
        vector<string> alignment = parseString(inp, "\t", 3);
        ++alignments;
        if (prevQName.compare(alignment[0])) {
            prevQName = alignment[0];
            ++reads;
        }
        if ((((uint) stoi(alignment[1])) & 0x4) == 0) {
            ++alignedAlignments;
        }
        if (chromosomeCounts.find(alignment[2]) == chromosomeCounts.end()) {
            chromosomeCounts.emplace(alignment[2], 1);
        } else {
            ++chromosomeCounts[alignment[2]];
        }
    }

    int mapped = 0;
    while (getline(in2, inp)) {
        mapped += stoi(parseString(inp, "\t", 3)[2]);
    }
   
    int time = -1;
    while(getline(in3, inp)) {
        vector<string> vec = parseString(inp, ":", 0);
        if (vec.size() != 3) {
            continue;
        }
        time = 3600 * stoi(vec[0]) + 60 * stoi(vec[1]) + stoi(vec[2]);
    }

    out << alignedAlignments << '\t' << alignments << '\t' << reads << '\t'
        << mapped << "\t" << chromosomeCounts.size() << '\t'
        << time << endl;

    in.close();
    in2.close();
    in3.close();
    out.close();
    return 0;
}

int categorizeUnmapped(string unmappedin, string output, bool sameQName,
        bool genomebam) {
    ifstream in(unmappedin);
    if (!in.is_open()) {
        cout << "Unable to open " << unmappedin << endl;
        return false;
    }
    ofstream out(output, ofstream::app);
    if (!out.is_open()) {
        cout << "Unable to open " << output << endl;
        return false;
    }
    string inp;
    unordered_set<string> qNames;
    while (getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '@') { continue; }
        vector<string> alignment = parseString(inp, "\t", 2);
        if (!sameQName) {
            alignment[0] = alignment[0].substr(0, alignment[0].size() - 2);
        }
        if (qNames.find(alignment[0]) != qNames.end()) { continue; }
        out << alignment[0] << "\t";
        uint flag = (uint) stoi(alignment[1]);
        bool unmapped = flag & 0x04;
        bool nextUnmapped = (flag & 0x01) && (flag & 0x08);
        if (unmapped || nextUnmapped) {
            if (genomebam) {
                if (unmapped && nextUnmapped) {
                    out << "unmapped" << endl;
                } else {
                    out << "mapped" << endl;
                }
            } else {
                out << "unaligned" << endl;
            }
        } else {
            out << "unmapped" << endl;
        }
        qNames.emplace(alignment[0]);
    }
    in.close();
    out.close();
    return true;
}

int getAllQNames(string sam, string outfile, bool sameQName) {
    ifstream in(sam);
    if (!in.is_open()) {
        cout << "Unable to open " << sam << endl;
        return false;
    }
    ofstream out(outfile);
    if (!out.is_open()) {
        cout << "Unable to open " << outfile << endl;
        return false;
    }
    string inp;
    unordered_set<string> allQNames;
    while(getline(in, inp)) {
        if (inp.size() == 0 || inp[0] == '@') { continue; }
        string qName = parseString(inp, "\t", 1)[0];
        if (!sameQName) {
            qName = qName.substr(0, qName.size() - 2);
        }
        if (allQNames.find(qName) == allQNames.end()) {
            out << qName << endl;
            allQNames.emplace(qName);
        }
    }
    in.close();
    out.close();
    return true;
}

int main(int argc, char **argv) {
    if (argc == 1) {
        cout << "no zeroes:    z infile outfile" << endl;
        cout << "count:        c infile" << endl;
        cout << "intersection: i infile1 infile2 outfile" << endl;
        cout << "stats:        s infile1 infile2" << endl;
        cout << "find:         f tofind transcript_num infile" << endl;
        cout << "ids:          k infile outfile transcriptome" << endl;
        cout << "get reads:    r insam ref_fastq outfastq" << endl;
        cout << "sparse:       x in_tsv out_tsv" << endl;
        cout << "names:        n infile outfile transcriptome" << endl;
        cout << "EQ to TCC:    t infile outprefix transcriptome ref_ec" << endl;
        cout << "t names:      w intranscriptome outtranscriptome" << endl;
        cout << "gene table:   g ingtf outtable" << endl;
        cout << "t index:      y transcriptome outtable" << endl;
        cout << "Check GFF:    q GFF flag" << endl;
        cout << "Timecourse:   b insam intsv inlog outtsv" << endl;
        cout << "Pull times:   u inlog" << endl;
        cout << "Unmapped cat  o unmappedin output sameQName genomebam" << endl;
        cout << "All QNAMEs    p insam output sameqName" << endl;
        return 1;
    }
    char opt = argv[1][1];
    
    int err = 1;
    switch (opt) {
        case 'z':   err = no_zeros(argv[2], argv[3]);
                    break;
        case 'c':   err = count(argv[2]);
                    break;
        case 'i':   err = intersection(argv[2], argv[3], argv[4]);
                    break;
        case 's':   err = stats(argv[2], argv[3]);
                    break;
        case 'f':   err = find_in_transcriptome_transcript(argv[2],
                                                           stoi(argv[3]),
                                                           argv[4]);
                    break;
        case 'k':   err = transcript_ids(argv[2], argv[3], argv[4]);
                    break;
        case 'r':   err = fill_fastq(argv[2], argv[3], argv[4]);
                    break;
        case 'x':   err = moar_zeroes(argv[2], argv[3]);
                    break;
        case 'n':   err = add_transcript_names(argv[2], argv[3], argv[4]);
                    break;
        case 't':   err = salmon_EQ_to_TCC(argv[2], argv[3], argv[4], argv[5]);
                    break;
        case 'w':   err = normal_transcriptome(argv[2], argv[3]);
                    break;
        case 'g':   err = gene_table(argv[2], argv[3]);
                    break;
        case 'y':   err = transcript_index(argv[2], argv[3]);
                    break;
        case 'q':   err = checkGFF(argv[2], (argc > 3) ? stoi(argv[3]) : 0,
                        (argc == 5) ? argv[4] : "");
                    break;
        case 'b':   err = timecourse(argv[2], argv[3], argv[4], argv[5]);
                    break;
        case 'u':   err = pull_times(argv[2]);
                    break;
        case 'o':   err = categorizeUnmapped(argv[2], argv[3],
                            argv[4][0] == '1', argv[5][0] == '1');
                    break;
        case 'p':   err = getAllQNames(argv[2], argv[3],
                            argv[4][0] == '1');
                    break;
    }
    
    return err;
}
