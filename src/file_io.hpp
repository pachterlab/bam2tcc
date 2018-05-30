/**
 * readGTFs and readSAM plus helper functions. See structs.hpp for more info
 * about the structs used.
 */

#ifndef file_io_util
#define file_io_util

#include <string>
#include <vector>

#include "structs.hpp"
#include "TCC_Matrix.hpp"

#define MAX_NOT_FOUND 25
#define MIN_UPDATE 1000000

void get_read(std::string info, Read &read);

Sequence get_sequence(std::string info);

int readGTFs(std::vector<std::string> files,
             std::vector<std::string> &transcriptome,
             std::vector<std::vector<Exon>*> &exons, int verbose);

int readSAM(std::string file, int filenumber,
            std::vector<std::vector<Exon>*> &exons, TCC_Matrix &matrix,
            std::string unmatched_outfile, int verbose);

#endif
