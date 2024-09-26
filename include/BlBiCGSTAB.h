#ifndef BLBICGSTAB_H
#define BLBICGSTAB_H

#include <LU.h>
#include <iostream>
#include <vector>
#include <algorithm>
#include <concepts>
#include <type_traits>
#include <fstream>
#include <chrono>
#include <fbinio.h>

// TODO:
// template<class MatrixT, class VectorT>
// concept  MatVecObj = requires(MatrixT mat, VectorT vec)
// {
//     requires std::is_same_v<typename MatrixT::Scalar, typename VectorT::Scalar>;    
//     mat+mat;
//     mat-mat;
//     mat*mat;
//     mat*vec;
//     vec+vec;
//     vec-vec;
//     vec*vec;
//     vec*mat;
//     mat.rows();
//     mat.cols();
//     vec.rows();
//     mat.matvec(vec);
//     mat.col(0);
//     mat.adjoint();
//     vec.adjoint();
//     LU_solve<MatrixT, VectorT>(mat, mat);
//     mat.trace();
    
// };

// template <class MatrixT>
// concept ScalarIsComplex = requires(MatrixT m) 
// {
// requires std::is_same_v<typename std::remove_cv_t<typename MatrixT::Scalar>, std::complex<typename std::remove_cv_t<typename MatrixT::Scalar>::value_type>>;
// }; 

template<class MatrixT, class VectorT>
class BlBiCGSTAB {
public:  
    MatrixT A;
    MatrixT B;
    MatrixT Xk;
    double epsilon;
    int N;
    int s;

    MatrixT Rk;
    MatrixT Pk;
    MatrixT R0c;
    MatrixT Vk;
    MatrixT alpha;
    MatrixT beta;
    MatrixT Sk;
    MatrixT Tk; 
    typename MatrixT::Scalar omega;

    double R0_2norm_max;

public:
    BlBiCGSTAB(MatrixT _A, MatrixT _B, MatrixT _X0,
               double epsilon);

    void solve();
    bool check_exit(MatrixT&);
    double norm2_max(MatrixT&);
};

#include <../src/BlBiCGSTAB.cpp>

#endif