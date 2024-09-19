#include <BiCGSTAB.h>
# include <fbinio.h>


BiCGSTAB::BiCGSTAB(Eigen::MatrixXd _A, Eigen::VectorXd _b, Eigen::VectorXd _x0, double _epsilon):
                    A(_A), b(_b), xk(_x0), epsilon(_epsilon)
{
    n = A.rows();
    vk = Eigen::VectorXd::Zero(n);
    alpha = 0;
    sk = Eigen::VectorXd::Zero(n);
    tk = Eigen::VectorXd::Zero(n);
    rkp1 = Eigen::VectorXd::Zero(n);
    omega = 0;
    beta = 0;
    rk = b - A*xk;
    r0c = rk;
    pk = rk;
}

void BiCGSTAB::solve(){
    int k = 0;
    int matvec_count = 0;
    std::ofstream logs("BiCGSTAB_logs.csv", std::ios::out | std::ios::trunc);
    logs << "k,time,res_2norm,matvec_count\n";
    auto start = std::chrono::high_resolution_clock::now();
    auto check = std::chrono::high_resolution_clock::now();
    while (k < n){
        vk = A*pk;
        matvec_count++;
        alpha = r0c.dot(rk)/r0c.dot(vk);
        sk = rk - alpha*vk;
        std::cout << "step: " << float(k) + 0.5 << ", (sk,sk)^1/2 = " << std::sqrt(sk.dot(sk)) << "\n\n";
        check = std::chrono::high_resolution_clock::now();
        logs << float(k) + 0.5 << ',' 
             << std::chrono::duration_cast<std::chrono::microseconds>(check-start).count() << "," 
             << std::sqrt(sk.dot(sk)) << ',' 
             << matvec_count << '\n';
        if (sk.dot(sk) < epsilon*epsilon){
            xk += alpha * pk;
            break;
        }
        tk = A*sk;
        matvec_count++;
        omega = tk.dot(sk)/tk.dot(tk);
        xk += alpha * pk + omega * sk;
        rkp1 = sk - omega * tk;
        std::cout << "step: " << k + 1 << ", (rk,rk)^1/2 = " << std::sqrt(rkp1.dot(rkp1)) << "\n\n" << "xk = \n" << xk << "\n\n\n\n";
        check = std::chrono::high_resolution_clock::now();
        logs << k + 1 << ',' 
             << std::chrono::duration_cast<std::chrono::microseconds>(check-start).count() << "," 
             << std::sqrt(rkp1.dot(rkp1)) << ',' 
             << matvec_count <<'\n';
        if (rkp1.dot(rkp1) < epsilon*epsilon)
            break;
        beta = r0c.dot(rkp1)/r0c.dot(rk)*alpha/omega;
        rk = rkp1;
        pk = rkp1 + beta * (pk - omega * vk);
        k+=1;
    }
    logs.close();
    write_binary("BiCGSTAB_solution.dat", xk);
    Eigen::MatrixXd sol(3, 3);
    read_binary("BiCGSTAB_solution.dat", sol);
    std::cout << sol << "\n\n";
}