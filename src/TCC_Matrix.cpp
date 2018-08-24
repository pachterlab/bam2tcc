#include "TCC_Matrix.hpp"
#include <fstream>
using namespace std;

/**
 * Constructer for new TCC_Matrix holding information for `file_count` number of
 * SAM files.
 * 
 * @param num_files    number of SAM files
 */
TCC_Matrix::TCC_Matrix(int file_count) {
    num_files = file_count;
    matrix = new unordered_map<string, int*>;
    sem = new Semaphore;
}

/**
 * Deconstructor for TCC_Matrix.
 */
TCC_Matrix::~TCC_Matrix() {
    for (auto it = matrix->begin(); it != matrix->end(); ++it) {
        delete[] it->second;
    }
    delete matrix;
    delete sem;
}

/**
 * Increments count for TCC in file number `file_num`.
 *
 * @param TCC         String representation of equivalence class of read.
 * @param file_num    Index of SAM file (should be less than num_files).
 */
void TCC_Matrix::inc_TCC(string TCC, int file_num) {
    sem->dec();
    auto it = matrix->find(TCC);
    if (it == matrix->end()) {
        int *arr = new int[num_files];
        for (int i = 0; i < num_files; ++i) {
            arr[i] = 0;
        }
        matrix->emplace(TCC, arr);
    }
    ++(*matrix)[TCC][file_num];
    sem->inc();
}

/**
 * Decrements count for TCC in file number <file_num>. Primarily for use in
 * the case that the SAM file multimaps the read (i.e. contains multiple entries
 * for it). Assumes TCC already exists in matrix.
 *
 * @precondition      TCC exists in matrix, i.e. has non-zero count for at least
 *                    one SAM file.
 * @throws            out_of_range error if TCC not found.
 * @param TCC         String representation of equivalence class of read.
 * @param file_num    Index of SAM file (should be less than num_files).
 */
void TCC_Matrix::dec_TCC(string TCC, int file_num) {
    sem->dec();
    --(*matrix)[TCC][file_num];
    sem->inc();
}

/**
 * Writes information in TCC_Matrix to files of names <outname>.ec and
 * <outname>.tsv. Returns 1 if error occurs in opening files, otherwise 0.
 *
 * @param outname    Name of output files (without file extension).
 * @return           1 if error occurs in opening files, otherwise 0.
 */
int TCC_Matrix::write_to_file(string outname, int num_transcripts) {
    ofstream ec(outname + ".ec");
    ofstream tsv(outname + ".tsv");
    if (!ec.is_open() || !tsv.is_open()) { return 1; }

    int count = 0;
    while (count < num_transcripts) {
        ec << count << '\t' << count << endl;
        tsv << count;
        auto it = matrix->find(to_string(count));
        if (it == matrix->end()) {
            for (int i = 0; i < num_files; ++i) {
                tsv << '\t' << 0;
            }
        } else {
            for (int i = 0; i < num_files; ++i) {
                tsv << '\t' << it->second[i];
            }
        }
        tsv << endl;
        ++count;
    }
    for (auto it = matrix->begin(); it != matrix->end(); ++it) {
        if (it->first.find(',') == string::npos
                && stoi(it->first) < num_transcripts) { continue; }
        ec << count << '\t' << it->first << endl;
        tsv << count;
        for (int i = 0; i < num_files; ++i) {
            tsv << '\t' << it->second[i];
        }
        tsv << endl;
        ++count;
    }
    ec.close();
    tsv.close();
    return 0;
}

/**
 * Writes information in TCC_Matrix to files of names <outname>.ec and
 * <outname>.tsv. Returns 1 if error occurs in opening files, otherwise 0.
 * Specifically, writes everything in kallisto sparse format.
 *
 * TODO: Is there a more cache-friendly way of doing this, or even implementing
 * this class in light of the fact that I have to output things this way?
 *
 * @param outname    Name of output files (without file extension).
 * @return           1 if error occurs in opening files, otherwise 0.
 */
int TCC_Matrix::write_to_file_sparse(string outname, int num_transcripts) {
    ofstream ec(outname + ".ec");
    ofstream tsv(outname + ".tsv");
    if (!ec.is_open() || !tsv.is_open()) { return 1; }

    int count = 0;
    while (count < num_transcripts) {
        ec << count << '\t' << count << endl;
        auto it = matrix->find(to_string(count));
        if (it != matrix->end()) {
            for (int j = 0; j < num_files; ++j) {
                if (it->second[j] != 0) {
                    tsv << count << '\t' << j << '\t' << it->second[j] << endl;
                }
            }
        }
        ++count;
    }
    for (auto it = matrix->begin(); it != matrix->end(); ++it) {
        if (it->first.find(',') == string::npos
                && stoi(it->first) < num_transcripts) { continue; }
        ec << count << '\t' << it->first << endl;
        for (int i = 0; i < num_files; ++i) {
            if (it->second[i] != 0) {
                tsv << count << '\t' << i << '\t' << it->second[i] << endl;
            }
        }
        ++count;
    }
    ec.close();
    tsv.close();
    return 0;
}

/**
 * Writes information in TCC_Matrix to files of names <outname>.ec and
 * <outname>.tsv in the order specified by vector <order>. Returns 1 if error
 * occurs in opening files, otherwise 0.
 *
 * @param outname       Name of output files (without file extension).
 *
 * @param order         Vector of strings describing the order in which to
 * output equivalence classes.
 *
 * @param ecs           Set of the strings in order. Used to speed up runtime.
 *
 * @return              1 if error occurs in opening files, otherwise 0.
 */
int TCC_Matrix::write_to_file_in_order(string outname,
                                       const vector<string> &order,
                                       const unordered_set<string> &ecs) {
    
    /* Open the file and die if something goes wrong */
    ofstream ec(outname + ".ec");
    ofstream tsv(outname + ".tsv");
    if (!ec.is_open() || !tsv.is_open()) {
        return 1;
    }
    
    /* Go through each equivalence class in vector order, use our map to locate
     the TCC, and output it */
    for (uint i = 0; i < order.size(); ++i) {
        ec << i << '\t' << order[i] << endl;
        
        // Try to find equivalence class in TCC matrix.
        auto elt = matrix->find(order[i]);
        tsv << i;
        // Equivalence class found. Add it to the vector of unfound classes.
        if (elt != matrix->end()) {
            for (int j = 0; j < num_files; ++j) {
                tsv << '\t' << elt->second[j];
            }
        }
        // Otherwise, just output 0s.
        else {
            for (int j = 0; j < num_files; ++j) {
                tsv << '\t' << 0;
            }
        }
        tsv << endl;
    }
    
    /* Write to file those classes that did not show up in kallisto */
    
    uint64_t index = order.size();
    // Go through our matrix and see if any of our classes wasn't in kallisto's
    // ec.
    for (auto it = matrix->begin(); it != matrix->end(); ++it) {
        auto elt = ecs.find(it->first);
        // It wasn't in kallisto's file, so it's not currently in our output.
        if (elt == ecs.end()) {
            ec << index << '\t' << it->first << endl;
            
            tsv << index;
            for (int j = 0; j < num_files; ++j) {
                tsv << '\t' << it->second[j];
            }
            tsv << endl;
            
            ++index;
        }
    }
    
    ec.close();
    tsv.close();
    return 0;
}


/**
 * Writes information in TCC_Matrix to files of names <outname>.ec and
 * <outname>.tsv in the order specified by vector <order>. Returns 1 if error
 * occurs in opening files, otherwise 0.
 *
 * @param outname       Name of output files (without file extension).
 *
 * @param order         Vector of strings describing the order in which to
 * output equivalence classes.
 *
 * @param ecs           Set of the strings in order. Used to speed up runtime.
 *
 * @return              1 if error occurs in opening files, otherwise 0.
 */
int TCC_Matrix::write_to_file_in_order_sparse(string outname,
                                       const vector<string> &order,
                                       const unordered_set<string> &ecs) {
    
    /* Open the file and die if something goes wrong */
    ofstream ec(outname + ".ec");
    ofstream tsv(outname + ".tsv");
    if (!ec.is_open() || !tsv.is_open()) {
        return 1;
    }

    /* Output all kallisto TCCs to ec file. */
    for (uint i = 0; i < order.size(); ++i) {
        ec << i << '\t' << order[i] << endl;
    }

    /* Go through each equivalence class in vector order, use our map to locate
     the TCC, and output it */
    for (int i = 0; i < num_files; ++i) {
        for (uint j = 0; j < order.size(); ++j) {
            // Try to find equivalence class in TCC matrix.
            auto elt = matrix->find(order[j]);
            // Equivalence class found and the count is not zero.
            if (elt != matrix->end() && elt->second[i] != 0) {
                tsv << j << '\t' << i << '\t' << elt->second[i];
                tsv << endl; 
            }
        }
    }
    
    /* Write to file those classes that did not show up in kallisto */
    uint64_t index = order.size();
    unordered_map<string, int> *output_index_map
                                            = new unordered_map<string, int>;
    // Go through our matrix and see if any of our classes weren't in kallisto's
    // ec.
    for (int i = 0; i < num_files; ++i) {
        for (auto it = matrix->begin(); it != matrix->end(); ++it) {
            auto elt = ecs.find(it->first);
            // It wasn't in kallisto's file so it's not currently in our output.
            if (elt == ecs.end() && it->second[i] != 0) {
                auto output_index =
                                    output_index_map->find(it->first);
                // If this TCC hasn't already been seen once, add it to our map
                // of TCC -> output index.
                if (output_index == output_index_map->end()) {
                    ec << index << '\t' << it->first << endl;
                    output_index_map->emplace(it->first, index);
                    ++index;
                }
                tsv << output_index_map->at(it->first) << '\t' << i << '\t';
                tsv << it->second[i] << endl;
            }
        }
    }
   
    delete output_index_map;
    ec.close(); 
    tsv.close();
    return 0;
}
