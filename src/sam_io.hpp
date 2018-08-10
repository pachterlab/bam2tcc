/**
 * Functions for reading SAM/BAM files.
 */
#ifndef __SAM_IO_HPP__
#define __SAM_IO_HPP__

#include <string>
#include <vector>
#include <unordered_map>
#include "exon.hpp"
#include "TCC_Matrix.hpp"
#include "semaphore.hpp"

#define PRIMARY_ONLY 0

int readSAM(std::string file, int filenumber,
             std::unordered_map<std::string, std::vector<Exon>*> &exons,
             TCC_Matrix &matrix,
             std::string unmatched_outfile, int verbose, int nthreads,
             bool rapmap, bool paired);

#endif
