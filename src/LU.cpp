#include <LU.h>

template<class MatrixT, class VectorT>
MatrixT LU_solve(MatrixT &A, MatrixT &B) 
{
    assert(A.cols()==A.rows());
    int s = B.cols();
    Eigen::FullPivLU<MatrixT> lu(A);
    MatrixT U = lu.matrixLU().template triangularView<Eigen::Upper>(); 
    int n = B.rows();
    MatrixT L =  MatrixT::Identity(n, n);
    L.block(0,0,n,n).template triangularView<Eigen::StrictlyLower>() = lu.matrixLU();
    MatrixT P = lu.permutationP();
    MatrixT Q = lu.permutationQ();
    // std::cout << L << "\n\n";
    // std::cout << U << "\n\n";
    // std::cout << P << "\n\n";
    // std::cout << Q << "\n\n";
    // std::cout << P.inverse()*L*U*Q.inverse() << "\n\n";
    MatrixT X(n, s);
    for (int i = 0; i < s; i++)
    {
        VectorT bw = P*B.col(i);
        VectorT xww = L_solve<MatrixT, VectorT>(L, bw);
        VectorT xw = U_solve<MatrixT, VectorT>(U, xww);
        X.col(i) = Q*xw;
    }
    return X;
}


template<class MatrixT, class VectorT>
VectorT L_solve(MatrixT &L, VectorT &b)
{
    typename VectorT::Scalar tmp;
    int n = L.rows();
    VectorT x(n);
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

template<class MatrixT, class VectorT>
VectorT U_solve(MatrixT &U, VectorT &b)
{
    typename VectorT::Scalar tmp;
    int n = U.rows();
    VectorT x(n);
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
