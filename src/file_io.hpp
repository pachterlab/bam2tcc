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
#include "semaphore.hpp"

#define NUM_GTF_ELT 9
#define NUM_SAM_ELT 11
#define ID_START "transcript_id \""
#define ID_END "\""
#define MAX_NOT_FOUND 25
#define MIN_UPDATE 1000000
#define USING_GTF 1

std::vector<int> get_eq(
        const std::vector<std::vector<Exon>*> &exons, Read &read);

void get_read(std::string info, Read &read);

Sequence get_sequence(std::string info);

int readSAMs(std::vector<std::string> &files,
             std::vector<std::vector<Exon>*> &exons, TCC_Matrix &matrix,
             std::string unmatched_outfile, int verbose, int nthreads);

int readGTFs(std::vector<std::string> &files,
             std::vector<std::string> &transcriptome,
             std::vector<std::vector<Exon>*> &exons, int verbose);

#endif
