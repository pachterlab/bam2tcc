#ifndef __COMMON_HPP__
#define __COMMON_HPP__

#include <string>
#include <vector>

std::vector<std::string> parseString(
        std::string query, std::string by, int to);

std::string lower(std::string s);

bool isNumber(std::string s);

#endif
