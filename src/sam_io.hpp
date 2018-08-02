/**
 * Functions for reading SAM/BAM files.
 */
#ifndef __SAM_IO_HPP__
#define __SAM_IO_HPP__

#include <string>
#include <vector>
#include "exon.hpp"
#include "TCC_Matrix.hpp"
#include "semaphore.hpp"

int readSAM(std::string file, int filenumber,
             std::vector<std::vector<Exon>*> &exons, TCC_Matrix &matrix,
             std::string unmatched_outfile, int verbose, int nthreads);

#endif
