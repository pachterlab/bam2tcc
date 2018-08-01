#ifndef __SAM_IO_HPP__
#define __SAM_IO_HPP__

#include <string>
#include <vector>

#include "exon.hpp"
#include "TCC_Matrix.hpp"
#include "semaphore.hpp"

#define MAX_NOT_FOUND 25
#define MIN_UPDATE 1000000
#define USING_GTF 1

int readSAM(std::string file, int filenumber,
             std::vector<std::vector<Exon>*> &exons, TCC_Matrix &matrix,
             std::string unmatched_outfile, int verbose, int nthreads);

#endif
