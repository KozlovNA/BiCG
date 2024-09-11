#include <BlBiCGSTAB.h>

BlBiCGSTAB::BlBiCGSTAB(Eigen::MatrixXd _A, Eigen::MatrixXd _B, Eigen::MatrixXd _X0, double _epsilon):
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
        R_norms[i] = Rk.col(i).transpose()*Rk.col(i);
    } 
    
    R0_2norm_max = *std::max_element(R_norms.begin(), R_norms.end());
}

void BlBiCGSTAB::solve()
{
    int k = 0;
    while (k < (N+s-1)/s){
        Vk = A*Pk;
        Eigen::MatrixXd alpha_system = R0c.transpose()*Vk;
        Eigen::MatrixXd alpha_rhs = R0c.transpose()*Rk;
        alpha = LU_solve(alpha_system, alpha_rhs);
        Sk = Rk - Vk*alpha;
        Tk = A*Sk;
        omega = (Tk.transpose()*Sk).trace() / ((Tk.transpose()*Tk).trace());
        Xk += Pk*alpha + omega*Sk;
        Rk = Sk - omega*Tk;
        if (check_exit())
        {
            break;
        }
        Eigen::MatrixXd beta_system = R0c.transpose()*Vk;
        Eigen::MatrixXd beta_rhs = -R0c.transpose()*Tk;
        beta = LU_solve(beta_system, beta_rhs);
        Eigen::MatrixXd Pk_prev = Pk;
        Pk = Rk + (Pk_prev - omega*Vk)*beta;
        k+=1;
    }
}

bool BlBiCGSTAB::check_exit()
{
    std::vector<double> R_norms(5, 0);
    for (int i = 0; i < s; i++) 
    {
        R_norms[i] = Rk.col(i).transpose()*Rk.col(i);
    } 
    double Rk_2norm_max = *std::max_element(R_norms.begin(), R_norms.end());
    return (Rk_2norm_max < (epsilon * epsilon * R0_2norm_max));
}