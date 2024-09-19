#ifndef FBINIO_H
#define FBINIO_H

# include <fstream>
# include <Eigen/Dense>

template<class Matrix>
void read_binary(const char* filename, Matrix& matrix);

template<class Matrix>
void write_binary(const char* filename, const Matrix& matrix);

# include "../src/fbinio.tpp"
#endif