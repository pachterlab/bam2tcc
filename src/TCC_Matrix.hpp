#ifndef tcc_matrix_hpp
#define tcc_matrix_hpp

#include <unordered_map>
#include <set>
#include <vector>
#include "semaphore.hpp"

/**
 * Class representing a matrix of TCC counts. Insertion and removal from the
 * matrix is thead-safe.
 */
class TCC_Matrix {
private:
    /* Number of SAM files data int this matrix represents */
    int num_files;
    /* Matrix holding data. Each int* is an array of size num_files */
    std::unordered_map<std::string, int*> *matrix;
    /* Semaphore to control access to this matrix. */
    Semaphore *sem;
public:
    TCC_Matrix(int num_files);
    ~TCC_Matrix();
    void inc_TCC(std::string TCC, int file_num);
    void dec_TCC(std::string TCC, int file_num);
    int write_to_file(std::string outname, int num_transcripts=0);
    int write_to_file_sparse(std::string outname);
    int write_to_file_in_order(std::string outname,
                               const std::vector<std::string> &order,
                               const std::set<std::string> &ecs);
    int write_to_file_in_order_sparse(std::string outname,
                               const std::vector<std::string> &order,
                               const std::set<std::string> &ecs);
};

#endif
