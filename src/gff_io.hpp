#ifndef GFF_IO_HPP
#define GFF_IO_HPP

#include <string>
#include <vector>
#include "structs.hpp"

#define NUM_GTF_ELT 9
#define ID_START "transcript_id \""
#define ID_END "\""
#define MAX_NOT_FOUND 25
#define MIN_UPDATE 1000000
#define USING_GTF 1

int readGFFs(std::vector<std::string> &files,
             std::vector<std::string> &transcriptome,
             std::vector<std::vector<Exon>*> &exons, int verbose);

#endif
