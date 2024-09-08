#ifndef BICGSTAB_H
#define BICGSTAB_H

#include<Eigen/Dense>
#include<iostream>

class BiCGSTAB {
public:
    int n;
    Eigen::MatrixXd A;
    Eigen::VectorXd b;
    Eigen::VectorXd r0c; 
    double epsilon; 

    Eigen::MatrixXd R;
    Eigen::MatrixXd V;
    Eigen::MatrixXd P;
    Eigen::VectorXd alpha;
    Eigen::MatrixXd S;
    Eigen::MatrixXd T;
    Eigen::VectorXd omega;
    Eigen::VectorXd beta;
    Eigen::VectorXd xk;
public:
    BiCGSTAB(Eigen::MatrixXd _A, Eigen::VectorXd _b, Eigen::VectorXd _x0,
             double _epsilon);

    void solve();
};

#endif