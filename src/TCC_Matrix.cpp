#include "TCC_Matrix.hpp"
#include <fstream>
using namespace std;

/**
 * Constructer for new TCC_Matrix holding information for <file_count> number of
 * SAM files.
 * 
 * @param num_files    number of SAM files
 */
TCC_Matrix::TCC_Matrix(int file_count) {
    num_files = file_count;
    matrix = new vector<int*>;
    indices = new unordered_map<string, int>;
}

/**
 * Deconstructor for TCC_Matrix.
 */
TCC_Matrix::~TCC_Matrix() {
    for (uint i = 0; i < matrix->size(); ++i) {
        delete[] (*matrix)[i];
    }
    delete matrix;
    delete indices;
}

/**
 * Increments count for TCC in file number <file_num>.
 *
 * @param TCC         String representation of equivalence class of read.
 * @param file_num    Index of SAM file (should be less than num_files).
 */
void TCC_Matrix::inc_TCC(string TCC, int file_num) {
    try {
        ++(*matrix)[indices->at(TCC)][file_num];
    }
    catch(out_of_range &e) {
        indices->emplace(TCC, matrix->size());
        int *arr = new int[num_files];
        for (int i = 0; i < num_files; ++i) {
            arr[i] = 0;
        }
        matrix->push_back(arr);
        ++(*matrix)[matrix->size() - 1][file_num];
    }
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
    --(*matrix)[indices->at(TCC)][file_num];
}

/**
 * Writes information in TCC_Matrix to files of names <outname>.ec and
 * <outname>.tsv. Returns 1 if error occurs in opening files, otherwise 0.
 *
 * @param outname    Name of output files (without file extension).
 * @return           1 if error occurs in opening files, otherwise 0.
 */
int TCC_Matrix::write_to_file(string outname) {
    ofstream ec(outname + ".ec");
    ofstream tsv(outname + ".tsv");
    if (!ec.is_open() || !tsv.is_open()) {
        return 1;
    }

    int count = 0;
    for (unordered_map<string, int>::iterator it = indices->begin();
        it != indices->end(); ++it) {

        ec << count << '\t' << it->first << endl;

        tsv << count;
        int *curr_TCC = (*matrix)[it->second];
        for (int i = 0; i < num_files; ++i) {
            tsv << '\t' << curr_TCC[i];
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
                                       const set<string> &ecs) {
    
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
        unordered_map<string, int>::iterator elt = indices->find(order[i]);
        tsv << i;
        // Equivalence class found. Add it to the vector of unfound classes.
        if (elt != indices->end()) {
            int *curr_TCC = (*matrix)[elt->second];
            for (int j = 0; j < num_files; ++j) {
                tsv << '\t' << curr_TCC[j];
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
    for (unordered_map<string, int>::iterator it = indices->begin();
         it != indices->end(); ++it) {
        
        set<string>::iterator elt = ecs.find(it->first);
        // It wasn't in kallisto's file, so it's not currently in our output.
        if (elt == ecs.end()) {
            ec << index << '\t' << it->first << endl;
            
            tsv << index;
            int *curr_TCC = (*matrix)[it->second];
            for (int j = 0; j < num_files; ++j) {
                tsv << '\t' << curr_TCC[j];
            }
            tsv << endl;
            
            ++index;
        }
    }
    
    tsv.close();
    return 0;
}
