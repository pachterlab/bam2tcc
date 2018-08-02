/**
 * Functions to generate kallisto-esque output, i.e. to map TCCs onto
 * kallisto's.
 */

#ifndef __KALLISTO_UTIL_HPP__
#define __KALLISTO_UTIL_HPP__

#include <vector>
#include <string>
#include <unordered_map>
#include <set>

//int get_index_to_seqid(const std::vector<std::string> &files,
//                       std::unordered_map<int, std::string> &map);

int get_index_to_kallisto_index(const std::vector<std::string> &gtf,
                            const std::vector<std::string> &transcriptome,
                            std::unordered_map<int, int> &map,
                            int verbose);

int change_index(const std::vector<std::string> &gtf,
               const std::vector<std::string> &transcriptome,
               std::string in_ec, std::string out_ec);

int get_kallisto_ec_order(std::string ec, std::vector<std::string> &v,
                          std::set<std::string> &s);

#endif
