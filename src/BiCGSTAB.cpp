#include <BiCGSTAB.h>

BiCGSTAB::BiCGSTAB(Eigen::MatrixXd _A, Eigen::VectorXd _b, Eigen::VectorXd _x0, double _epsilon):
                    A(_A), b(_b), xk(_x0), epsilon(_epsilon)
{
    n = A.rows();
    vk = Eigen::VectorXd::Zero(n);
    alpha = 0;
    sk = Eigen::VectorXd::Zero(n);
    tk = Eigen::VectorXd::Zero(n);
    omega = 0;
    beta = 0;
    rk = b - A*xk;
    r0c = rk;
    pk = rk;
}

void BiCGSTAB::solve(){
    int k = 0;
    while (k < n){
        vk = A*pk;
        alpha = r0c.dot(rk)/r0c.dot(vk);
        sk = rk - alpha*vk;
        std::cout << "step: " << k << ", (sk,sk)^1/2 = " << std::sqrt(sk.dot(sk)) << "\n\n";
        tk = A*sk;
        omega = tk.dot(sk)/tk.dot(tk);
        xk += alpha * pk + omega * sk;
        if (sk.dot(sk) < epsilon*epsilon)
            break;
        Eigen::VectorXd rkp1 = sk - omega * tk;
        beta = r0c.dot(rkp1)/r0c.dot(rk)*alpha/omega;
        rk = rkp1;
        pk = rkp1 + beta * (pk - omega * vk);
        k+=1;
    }
    std::cout << "step: " << k << ", (sk,sk)^1/2 = " << std::sqrt(sk.dot(sk)) << "\n\n";
}