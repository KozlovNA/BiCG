#ifndef LU_H
#define LU_H
#include <Eigen/Dense>
#include <iostream>

template<class MatrixT, class VectorT>
MatrixT LU_solve(MatrixT &A, MatrixT &B); 
    // Solves linear system with many right-hand sides AX=B  
template<class MatrixT, class VectorT>          
VectorT L_solve(MatrixT &L, VectorT &b);

template<class MatrixT, class VectorT>
VectorT U_solve(MatrixT &U, VectorT &b);

#include <../src/LU.cpp>
#endif