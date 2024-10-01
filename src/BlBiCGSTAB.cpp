#include <BlBiCGSTAB.h>

template<class MatrixT, class VectorT>
BlBiCGSTAB<MatrixT, VectorT>::BlBiCGSTAB(MatrixT  _A, MatrixT  _B, MatrixT  _X0, double _epsilon):
                    A(_A), B(_B), Xk(_X0), epsilon(_epsilon)
{
    N = Xk.rows();
    s = Xk.cols();
    Rk = B - A*Xk;
    R0c = Rk;
    // Pk = Rk;
    Pk = orth(Rk, epsilon);
    
    R0_2norm_max = norm2_max(R0c);
}

template<class MatrixT, class VectorT>
void BlBiCGSTAB<MatrixT, VectorT>::solve()
{
    int k = 0;
    int matvec_count = 0;
    std::ofstream logs("../output/BlBiCGSTAB_logs.csv", std::ios::out | std::ios::trunc);
    logs << "k,time,res_max2norm_rel,matvec_count\n";
    auto start = std::chrono::high_resolution_clock::now();
    auto check = std::chrono::high_resolution_clock::now();

    while (k < (N+s-1)/s){
        Vk = A.matvec(Pk);
        matvec_count += Pk.cols();
        MatrixT alpha_system = R0c.adjoint()*Vk;
        MatrixT alpha_rhs = R0c.adjoint()*Rk;
        alpha = QR_solve<MatrixT, VectorT>(alpha_system, alpha_rhs);
        Sk = Rk - Vk*alpha;

        double sk_norm_sq_rel = norm2_max(Sk)/R0_2norm_max;      
        std::cout << "step: " << float(k) + 0.5 << ", residuals 2norm ratio = " << std::sqrt(sk_norm_sq_rel) << "\n\n";
        check = std::chrono::high_resolution_clock::now();
        logs << float(k) + 0.5 << ',' 
             << std::chrono::duration_cast<std::chrono::microseconds>(check-start).count() << "," 
             << std::sqrt(sk_norm_sq_rel) << ',' 
             << matvec_count << '\n';
        
        if (check_exit(Sk)){
           Xk += Pk*alpha;
           break;}

        Tk = A.matvec(Sk);
        matvec_count += Sk.cols();
        omega = (Tk.adjoint()*Sk).trace() / ((Tk.adjoint()*Tk).trace());
        Xk += Pk*alpha + omega*Sk;
        Rk = Sk - omega*Tk;

        double rk_norm_sq_rel = norm2_max(Rk)/R0_2norm_max;      
        std::cout << "step: " << float(k) + 1 << ", residuals 2norm ratio = " << std::sqrt(rk_norm_sq_rel) << "\n\n";
        check = std::chrono::high_resolution_clock::now();
        logs << float(k) + 1 << ',' 
             << std::chrono::duration_cast<std::chrono::microseconds>(check-start).count() << "," 
             << std::sqrt(rk_norm_sq_rel) << ',' 
             << matvec_count << '\n';

        if (check_exit(Rk)) break;

        MatrixT beta_system = R0c.adjoint()*Vk;
        MatrixT beta_rhs = -R0c.adjoint()*Tk;
        beta = QR_solve<MatrixT, VectorT>(beta_system, beta_rhs);
        MatrixT Pk_prev = Pk;
        MatrixT Pkt = Rk + (Pk_prev - omega*Vk)*beta;
        Pk = orth(Pkt, epsilon);
        k+=1;
    }
    logs.close();
    write_binary("../output/BlBiCGSTAB_solution.dat", Xk);
    std::cout << Xk << "\n\n";
}

template<class MatrixT, class VectorT>
bool BlBiCGSTAB<MatrixT, VectorT>::check_exit(MatrixT& Res)
{
    double Res_2norm_max = norm2_max(Res);
    return (Res_2norm_max < (epsilon * epsilon * R0_2norm_max));
}

template<class MatrixT, class VectorT>
double BlBiCGSTAB<MatrixT, VectorT>::norm2_max(MatrixT& Mat)
{
    int n = Mat.cols();
    std::vector<double> R_norms(n, 0);
    for (int i = 0; i < n; i++) 
    {
        R_norms[i] = Mat.col(i).squaredNorm();
    } 
    double Res_2norm_max = *std::max_element(R_norms.begin(), R_norms.end());
    return Res_2norm_max;
}

template<class MatrixT, class VectorT>
MatrixT BlBiCGSTAB<MatrixT, VectorT>::orth(MatrixT& A, double eps)
{
    MatrixT Q;
    // MatrixT P;
    Eigen::ColPivHouseholderQR<Eigen::MatrixX<typename MatrixT::Scalar>> qr(A);
    // std::cout << qr.rank() << "\n\n";
    qr.setThreshold(eps);
    Q = qr.householderQ();
    // P = qr.colsPermutation();
    return Q.leftCols(qr.rank());
}

template<class MatrixT, class VectorT>
MatrixT QR_solve(MatrixT& A, MatrixT& B)
{
    MatrixT Q;
    MatrixT R;
    Eigen::HouseholderQR<Eigen::MatrixX<typename MatrixT::Scalar>> qr(A);
    Q = qr.householderQ();
    R = qr.matrixQR().template triangularView<Eigen::Upper>();
    MatrixT B1 = Q.adjoint()*B;
    MatrixT X(A.rows(), B.cols());
    for (int i = 0; i < B1.cols(); i++)
    {
        VectorT b = B1.col(i);
        X.col(i) = U_solve(R, b);
    }  
    return X;
}