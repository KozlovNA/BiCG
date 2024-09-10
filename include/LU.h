#ifndef LU_H
#define LU_H
#include <Eigen/Dense>
#include <iostream>

Eigen::MatrixXd LU_solve(Eigen::MatrixXd &A, Eigen::MatrixXd &B); 
    // Solves linear system with many right-hand sides AX=B        
Eigen::VectorXd L_solve(Eigen::MatrixXd &L, Eigen::VectorXd &b);
Eigen::VectorXd U_solve(Eigen::MatrixXd &U, Eigen::VectorXd &b);

#endif