#ifndef MATVEC_H
#define MATVEC

#include <Eigen/Dense>

class MyMatrixXcf: public Eigen::MatrixXcf 
{
public:
    using Eigen::MatrixXcf::MatrixXcf;
    Eigen::MatrixXcf matvec(Eigen::MatrixXcf& B) 
    {
        return (*this)*B;
    }
};

class MyMatrixXd: public Eigen::MatrixXd 
{
public:
    using Eigen::MatrixXd::MatrixXd;
    Eigen::MatrixXd matvec(Eigen::MatrixXd& B) 
    {
        return (*this)*B;
    }
};

#endif