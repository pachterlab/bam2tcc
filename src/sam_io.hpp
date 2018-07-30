#ifndef SAM_IO_HPP
#define SAM_IO_HPP

#include <string>
#include <vector>

#include "structs.hpp"
#include "TCC_Matrix.hpp"
#include "semaphore.hpp"

#define MAX_NOT_FOUND 25
#define MIN_UPDATE 1000000
#define USING_GTF 1

std::vector<int> get_eq(
        const std::vector<std::vector<Exon>*> &exons, Read &read);

void get_read(std::string info, Read &read);

int readSAMs(std::vector<std::string> &files,
             std::vector<std::vector<Exon>*> &exons, TCC_Matrix &matrix,
             std::string unmatched_outfile, int verbose, int nthreads);

#endif
