#ifndef FBINIO_H
#define FBINIO_H

# include <fstream>
# include <Eigen/Dense>

template<class Matrix>
void read_binary(const char* filename, Matrix& matrix){
    std::ifstream in(filename, std::ios::in | std::ios::binary);
    uint32_t rows=0, cols=0;
    in.read((char*) (&rows),sizeof(uint32_t));
    in.read((char*) (&cols),sizeof(uint32_t));
    matrix.resize(rows, cols);
    in.read((char *) matrix.data() , rows*cols*sizeof(typename Matrix::Scalar));
    in.close();
}
template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix){
    std::ofstream out(filename, std::ios::out | std::ios::binary | std::ios::trunc);
    uint32_t rows=matrix.rows(), cols=matrix.cols();
    out.write((char*) (&rows), sizeof(uint32_t));
    out.write((char*) (&cols), sizeof(uint32_t));
    out.write((char*) matrix.data(), rows*cols*sizeof(typename Matrix::Scalar) );
    out.close();
}
#endif