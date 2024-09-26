#include <BlBiCGSTAB.h>

template<class MatrixT, class VectorT>
BlBiCGSTAB<MatrixT, VectorT>::BlBiCGSTAB(MatrixT  _A, MatrixT  _B, MatrixT  _X0, double _epsilon):
                    A(_A), B(_B), Xk(_X0), epsilon(_epsilon)
{
    N = Xk.rows();
    s = Xk.cols();
    Rk = B - A*Xk;
    Pk = Rk;
    R0c = Rk;

    std::vector<double> R_norms(s, 0);
    for (int i = 0; i < s; i++) 
    {
        R_norms[i] = Rk.col(i).squaredNorm();
    } 
    
    R0_2norm_max = *std::max_element(R_norms.begin(), R_norms.end());
}

template<class MatrixT, class VectorT>
void BlBiCGSTAB<MatrixT, VectorT>::solve()
{
    int k = 0;
    while (k < (N+s-1)/s){
        Vk = A*Pk;
        MatrixT alpha_system = R0c.adjoint()*Vk;
        MatrixT alpha_rhs = R0c.adjoint()*Rk;
        alpha = LU_solve<MatrixT, VectorT>(alpha_system, alpha_rhs);
        Sk = Rk - Vk*alpha;
        Tk = A*Sk;
        omega = (Tk.adjoint()*Sk).trace() / ((Tk.adjoint()*Tk).trace());
        Xk += Pk*alpha + omega*Sk;
        Rk = Sk - omega*Tk;
        if (check_exit()) break;
        MatrixT beta_system = R0c.adjoint()*Vk;
        MatrixT beta_rhs = -R0c.adjoint()*Tk;
        beta = LU_solve<MatrixT, VectorT>(beta_system, beta_rhs);
        MatrixT Pk_prev = Pk;
        Pk = Rk + (Pk_prev - omega*Vk)*beta;
        k+=1;
    }
}

// template<class VectorT>
// requires ScalarIsComplex<VectorT>
// double norm_sq(VectorT v)
// {
//     return (v.dot(v)).real();
// } 
// template <class VectorT>
// double norm_sq(VectorT v)
// {
//     return v.dot(v);
// } 

template<class MatrixT, class VectorT>
bool BlBiCGSTAB<MatrixT, VectorT>::check_exit()
{
    std::vector<double> R_norms(s, 0);
    for (int i = 0; i < s; i++) 
    {
        R_norms[i] = Rk.col(i).norm();
    } 
    double Rk_2norm_max = *std::max_element(R_norms.begin(), R_norms.end());
    return (Rk_2norm_max < (epsilon * epsilon * R0_2norm_max));
}