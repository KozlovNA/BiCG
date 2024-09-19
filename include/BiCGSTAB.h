#ifndef BICGSTAB_H
#define BICGSTAB_H


#include <Eigen/Dense>
#include <iostream>
#include <fstream>
#include <chrono>

class BiCGSTAB {
public:
    int n;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd r0c; 
    double epsilon; 

    Eigen::VectorXd rk;
    Eigen::VectorXd vk;
    Eigen::VectorXd pk;
    Eigen::VectorXd rkp1;
    double alpha;
    Eigen::VectorXd sk;
    Eigen::VectorXd tk;
    double omega;
    double beta;
    Eigen::VectorXd xk;
public:
    BiCGSTAB(Eigen::MatrixXd _A, Eigen::VectorXd _b, Eigen::VectorXd _x0,
             double _epsilon);

    void solve();

};

#endif