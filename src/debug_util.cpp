#include <fstream>
#include <iostream>
#include <unordered_map>

#include "util.hpp"
using namespace std;

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
        vector<string> tcc = parse_tsv(inptsv);
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
        vector<string> tsv = parse_tsv(inp);
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
        cerr << "ERROR: unable to open" << inname << ".tsv" << endl;
        return 1;
    }
    
    string inp;
    getline(in, inp);
    vector<string> tcc = parse_tsv(inp);
    
    vector<uint> count;
    for (uint i = 1; i < tcc.size(); ++i) {
        count.push_back((uint) stoi(tcc[i]));
    }
    
    while (getline(in, inp)) {
        tcc = parse_tsv(inp);
        for (uint i = 1; i < tcc.size(); ++i) {
            count[i - 1] += (uint) (stoi(tcc[i]));
        }
    }
    
    for (uint i = 0; i < tcc.size() - 1; ++i) {
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
    int index1 = stoi(parse_tsv(inp1)[0]);
    int index2 = stoi(parse_tsv(inp2)[0]);
    while (true) {
        if (index1 == index2) {
            out << inp1 << endl << inp2 << endl << endl;
            if (!getline(in1, inp1)) {
                break;
            }
            index1 = stoi(parse_tsv(inp1)[0]);
            if (!getline(in2, inp2)) {
                break;
            }
            index2 = stoi(parse_tsv(inp2)[0]);
        }
        else if (index1 < index2) {
            if (!getline(in1, inp1)) {
                break;
            }
            index1 = stoi(parse_tsv(inp1)[0]);
        }
        else {
            if (!getline(in2, inp2)) {
                break;
            }
            index2 = stoi(parse_tsv(inp2)[0]);
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
        int count1 = stoi(parse_tsv(inp1)[1]);
        int count2 = stoi(parse_tsv(inp2)[1]);
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
        uint first_space = inp.find('.');
        // Exclude the first character ('>') and the END CHAR.
        map->emplace(to_string(count), inp.substr(1, first_space - 1));
        ++count;
    }
    
    while (getline(in, inp)) {
        inp = lower(inp);
        vector<string> t1 = parse_tsv(inp);
        out << t1[0] << '\t';
        vector<string> t2 = parse_csv(t1[1]);
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

int main(int argc, char **argv) {
    if (argc == 1) {
        cout << "no zeroes:    z infile outfile" << endl;
        cout << "count:        c infile" << endl;
        cout << "intersection: i infile1 infile2 outfile" << endl;
        cout << "stats:        s infile1 infile2" << endl;
        cout << "find:         f tofind transcript_num infile" << endl;
        cout << "ids:          k infile outfile transcriptome" << endl;
        cout << "get reads     r insam ref_fastq outfastq" << endl;
        cout << "sparse        x in_tsv out_tsv" << endl;
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
    }
    
    return err;
}
