#ifndef BLBICGSTAB_H
#define BLBICGSTAB_H

#include<Eigen/Dense>
#include<LU.h>
#include<iostream>

class BlBiCGSTAB {
public:
    Eigen::MatrixXd A;
    Eigen::MatrixXd B;
    Eigen::MatrixXd Xk;
    double epsilon;
    int N;
    int s;

    Eigen::MatrixXd Rk;
    Eigen::MatrixXd Pk;
    Eigen::MatrixXd R0c;
    Eigen::MatrixXd Vk;
    Eigen::MatrixXd alpha;
    Eigen::MatrixXd beta;
    Eigen::MatrixXd Sk;
    Eigen::MatrixXd Tk; 
    double omega;


public:
    BlBiCGSTAB(Eigen::MatrixXd _A, Eigen::MatrixXd _B, Eigen::MatrixXd _X0,
               double epsilon);

    void solve();
};

#endif