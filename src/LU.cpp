#include <LU.h>

Eigen::MatrixXd LU_solve(Eigen::MatrixXd &A, Eigen::MatrixXd &B) 
{
    assert(A.cols()==A.rows());
    int s = B.cols();
    Eigen::FullPivLU<Eigen::MatrixXd> lu(A);
    Eigen::MatrixXd U = lu.matrixLU().triangularView<Eigen::Upper>(); 
    int n = B.rows();
    Eigen::MatrixXd L =  Eigen::MatrixXd::Identity(n, n);
    L.block(0,0,n,n).triangularView<Eigen::StrictlyLower>() = lu.matrixLU();
    Eigen::MatrixXd P = lu.permutationP();
    Eigen::MatrixXd Q = lu.permutationQ();
    // std::cout << L << "\n\n";
    // std::cout << U << "\n\n";
    // std::cout << P << "\n\n";
    // std::cout << Q << "\n\n";
    // std::cout << P.inverse()*L*U*Q.inverse() << "\n\n";
    Eigen::MatrixXd X(n, s);
    for (int i = 0; i < s; i++)
    {
        Eigen::VectorXd bw = P*B.col(i);
        Eigen::VectorXd xww = L_solve(L, bw);
        Eigen::VectorXd xw = U_solve(U, xww);
        X.col(i) = Q*xw;
    }
    return X;
}



Eigen::VectorXd L_solve(Eigen::MatrixXd &L, Eigen::VectorXd &b)
{
    double tmp;
    int n = L.rows();
    Eigen::VectorXd x(n);
    for (int i = 0; i < n; i++)
    {
        tmp = b(i);
        for (int j = 0; j < i; j++)
        {
            tmp -= L(i,j) * x(j);
        }
        x(i) = tmp / L(i, i);
    }
    return x;
}

Eigen::VectorXd U_solve(Eigen::MatrixXd &U, Eigen::VectorXd &b)
{
    double tmp;
    int n = U.rows();
    Eigen::VectorXd x(n);
    for (int i = n-1; i > -1; i--)
    {
        tmp = b(i);
        for (int j = i+1; j < n; j++)
        {
            tmp -= U(i,j) * x(j);
        }
        x(i) = tmp / U(i, i);
    }
    return x;
}
