#include "common.hpp"
using namespace std;

vector<string> parseString(string query, string by, int to) {
    vector<string> parsed;
    int start = 0, end = query.find(by);
    while (end != string::npos) {
        parsed.push_back(query.substr(start, end - start));
        if (to != 0 && parsed.size() == to) { return parsed; }
        start = end + by.length();
        end = query.find(by, start);
    }
    parsed.push_back(query.substr(start, query.length() - end));
    return parsed;
}

/**
 * @brief Returns copy of s with all alphabetical characters in lowercase.
 *
 * @param s     string to lower
 *
 * @return      copy of s in lowercase
 */
string lower(string s) {
    string l = "";
    for (int i = 0; i < s.length(); ++i)
        l += tolower(s[i]);
    return l;
}

bool isNumber(string s) {
    return s.find_first_not_of("0123456789") == string::npos;
}
