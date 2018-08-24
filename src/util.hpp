/**
 * Various general utility functions.
 */

#ifndef __UTIL_HPP__
#define __UTIL_HPP__

#include <string>
#include <vector>

std::string lower(std::string s);

uint64_t get_line_count(std::string filename);

int is_number(std::string s);

std::vector<std::string> parse_by(std::string query, std::string by);

std::vector<std::string> parse_tsv(std::string tsv);

std::vector<std::string> parse_csv(std::string csv);

std::vector<uint64_t> parse_csv_ints(std::string csv);

int test_open(std::string filename, int mode = 0);

#endif

