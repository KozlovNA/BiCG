#include <BiCGSTAB.h>

BiCGSTAB::BiCGSTAB(Eigen::MatrixXd _A, Eigen::VectorXd _b, Eigen::VectorXd _x0, double _epsilon):
                    A(_A), b(_b), xk(_x0), epsilon(_epsilon)
{
    n = A.rows();
    R = Eigen::MatrixXd::Zero(n, n+1);
    V = Eigen::MatrixXd::Zero(n, n+1);
    P = Eigen::MatrixXd::Zero(n, n+1);
    alpha = Eigen::VectorXd::Zero(n);
    S = Eigen::MatrixXd::Zero(n, n+1);
    T = Eigen::MatrixXd::Zero(n, n+1);
    omega = Eigen::VectorXd::Zero(n+1);
    beta = Eigen::VectorXd::Zero(n+1);
    R.col(0) = b - A*xk;
    r0c = R.col(0);
    P.col(0) = R.col(0);
}

void BiCGSTAB::solve(){
    int k = 0;
    while (k < n){
        V.col(k) = A*P.col(k);
        alpha(k) =  r0c.dot(R.col(k))/r0c.dot(V.col(k));
        S.col(k) = R.col(k) - alpha(k)*V.col(k);
        std::cout << "step: " << k << ", (sk,sk)^1/2 = " << std::sqrt(S.col(k).dot(S.col(k))) << "\n\n";
        T.col(k) = A*S.col(k);
        omega(k) = T.col(k).dot(S.col(k))/T.col(k).dot(T.col(k));
        xk += alpha(k) * P.col(k) + omega(k) * S.col(k);
        if (S.col(k).dot(S.col(k)) < epsilon*epsilon)
            break;
        R.col(k+1) = S.col(k) - omega(k) * T.col(k);
        beta(k) = r0c.dot(R.col(k+1))/r0c.dot(R.col(k))*alpha(k)/omega(k);
        P.col(k+1) = R.col(k+1) + beta(k) * (P.col(k) - omega(k) * V.col(k));
        k+=1;
    }
    std::cout << "step: " << k << ", (sk,sk)^1/2 = " << std::sqrt(S.col(k).dot(S.col(k))) << "\n\n";
}