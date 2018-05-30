#ifndef tcc_matrix_hpp
#define tcc_matrix_hpp

#include <unordered_map>
#include <set>
#include <vector>

/**
 * Class representing a matrix of TCC counts.
 */
class TCC_Matrix {
private:
    /* number of SAM files data int this matrix represents */
    int num_files;
    /* matrix holding data. each int * is an array of size num_files */
    std::vector<int*> *matrix;
    /* maps equivalence class in string form to index of EC in matrix */
    std::unordered_map<std::string, int> *indices;
public:
    TCC_Matrix(int num_files);
    ~TCC_Matrix();
    void inc_TCC(std::string TCC, int file_num);
    void dec_TCC(std::string TCC, int file_num);
    int write_to_file(std::string outname);
    int write_to_file_in_order(std::string outname,
                               const std::vector<std::string> &order,
                               const std::set<std::string> &ecs);
};

#endif
