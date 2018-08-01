/**
 * Various general utility functions.
 */

#include "util.hpp"
#include <fstream>      // file i/o
using namespace std;

/**
 * @brief Returns copy of s with all alphabetical characters in lowercase.
 *
 * @param s     string to lower
 *
 * @return      copy of s in lowercase
 */
string lower(string s) {
    string l = "";
    for (uint i = 0; i < s.length(); ++i)
        l += tolower(s[i]);
    return l;
}

/**
 * @brief Get line count of file. Since we only increment a counter, it is
 * reasonably fast.
 *
 * @param filename      name of file
 */
uint64_t get_line_count(string filename) {
    string inp;
    uint64_t line_count = 0;
    ifstream f(filename);
    if (!f.is_open()) {
        return -1;
    }
    
    while (getline(f, inp)) {
        ++line_count;
    }
    
    return line_count;
}

/**
 * @brief Tests whether a string is a number. Returns true if it contains only
 * digits.
 *
 * @param s     string to test
 * @return      true if string represents a number, else false
 */
int is_number(string s) {
    return s.find_first_not_of("0123456789") == string::npos;
}

/**
 * @brief Parses a tab-separated string and returns fields as elements in a
 * vector.
 *
 * @param tsv       string to parse
 * @return          vector containing tab-separated fields of string
 */
vector<string> parse_tsv(string tsv) {
    vector<string> values;
    int start = 0, end = tsv.find('\t');
    while (end != string::npos) {
        values.push_back(tsv.substr(start, end - start));
        start = end + 1;
        end = tsv.find('\t', start);
    }
    values.push_back(tsv.substr(start, tsv.length() - end));
    return values;
}

/**
 * @brief Parses a comma-separated string and returns fields as elements in a
 * vector.
 *
 * @param csv       string to parse
 * @return          vector containing comma-separated fields of string
 */
vector<string> parse_csv(string csv) {
    vector<string> values;
    int start = 0, end = csv.find(',');
    while (end != string::npos) {
        values.push_back(csv.substr(start, end - start));
        start = end + 1;
        end = csv.find(',', start);
    }
    values.push_back(csv.substr(start, csv.length() - end));
    return values;
}

/**
 * @brief Parses a comma-separated string and returns fields as elements in a
 * vector. This one assumes that all elements are ints, and will return a
 * vector of uint64_t types.
 *
 * @param csv       string to parse
 * @return          vector containing comma-separated fields of string
 */
vector<uint64_t> parse_csv_ints(string csv) {
    vector<uint64_t> values;
    int start = 0, end = csv.find(',');
    while (end != string::npos) {
        values.push_back(stoi(csv.substr(start, end - start)));
        start = end + 1;
        end = csv.find(',', start);
    }
    values.push_back(stoi(csv.substr(start, csv.length() - end)));
    return values;
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
int test_open(string filename, int mode) {
    if (mode == 0) {
        ifstream f(filename);
        if (f.is_open()) {
            f.close();
            return 1;
        }
        else {
            return 0;
        }
    }
    else if (mode == 1) {
        ofstream f(filename);
        if (f.is_open()) {
            f.close();
            return 1;
        }
        else {
            return 0;
        }
    }
    else if (mode == 2) {
        fstream f(filename);
        if (f.is_open()) {
            f.close();
            return 1;
        }
        else {
            return 0;
        }
    }
    else {
        return 0;
    }
}
