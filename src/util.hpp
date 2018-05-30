/**
 * Various general utility functions.
 */

#ifndef util_hpp
#define util_hpp

#include <string>
#include <vector>

std::string lower(std::string s);

uint64_t get_line_count(std::string filename);

int is_number(std::string s);

std::vector<std::string> parse_tsv(std::string tsv);

std::vector<std::string> parse_csv(std::string csv);

std::vector<uint64_t> parse_csv_ints(std::string csv);

int test_open(std::string filename, int mode = 0);

#endif /* common_h */
