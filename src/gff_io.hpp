/**
 * Functions for reading GFF files.
 */
#ifndef __GFF_IO_HPP__
#define __GFF_IO_HPP__

#include <string>
#include <vector>
#include <unordered_map>
#include "exon.hpp"

int readGFFs(std::vector<std::string> &files,
             std::vector<std::string> &transcriptome,
             std::unordered_map<std::string, std::vector<Exon>*> &exons,
             int verbose);

#endif
