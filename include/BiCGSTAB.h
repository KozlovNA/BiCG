#ifndef BICGSTAB_H
#define BICGSTAB_H


#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <chrono>
#include <fbinio.h>

template<typename T>
class BiCGSTAB {
public:
    typedef Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> MatrixXT;
    typedef Eigen::Matrix<T, Eigen::Dynamic, 1> VectorXT;
    int n;
    MatrixXT A;
    VectorXT b;
    VectorXT r0c; 
    double epsilon; 

    VectorXT rk;
    VectorXT vk;
    VectorXT pk;
    VectorXT rkp1;
    T alpha;
    VectorXT sk;
    VectorXT tk;
    T omega;
    T beta;
    VectorXT xk;
public:
    BiCGSTAB(MatrixXT _A, VectorXT _b, VectorXT _x0,
             double _epsilon);

    void solve();

    double norm_sq (VectorXT vec);



};

# include <../src/BiCGSTAB.cpp>

#endif