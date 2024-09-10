#include <BlBiCGSTAB.h>

BlBiCGSTAB::BlBiCGSTAB(Eigen::MatrixXd _A, Eigen::MatrixXd _B, Eigen::MatrixXd _X0, double _epsilon):
                    A(_A), B(_B), Xk(_X0), epsilon(_epsilon)
{
    N = Xk.rows();
    s = Xk.cols();
    Rk = B - A*Xk;
    Pk = Rk;
    R0c = Rk;
}

void BlBiCGSTAB::solve(){
    int k = 0;
    while (k < N*s){
        Vk = A*Pk;
        Eigen::MatrixXd alpha_system = R0c.transpose()*Vk;
        Eigen::MatrixXd alpha_rhs = R0c.transpose()*Rk;
        alpha = LU_solve(alpha_system, alpha_rhs);
        Sk = Rk - Vk*alpha;
        Tk = A*Sk;
        omega = (Tk.transpose()*Sk).trace() / ((Tk.transpose()*Tk).trace());
        Xk += Pk*alpha + omega*Sk;
        Rk = Sk - omega*Tk;
        Eigen::MatrixXd beta_system = R0c.transpose()*Vk;
        Eigen::MatrixXd beta_rhs = -R0c.transpose()*Tk;
        beta = LU_solve(beta_system, beta_rhs);
        Eigen::MatrixXd Pk_prev = Pk;
        Pk = Rk + (Pk_prev - omega*Vk)*beta;
        k+=1;
    }
}