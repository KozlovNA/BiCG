#include <BiCGSTAB.h>

template<typename T>
BiCGSTAB<T>::BiCGSTAB(MatrixXT _A, VectorXT _b, VectorXT _x0, double _epsilon):
                    A(_A), b(_b), xk(_x0), epsilon(_epsilon)
{
    n = A.rows();
    vk = VectorXT::Zero(n);
    alpha = 0;
    sk = VectorXT::Zero(n);
    tk = VectorXT::Zero(n);
    rkp1 =VectorXT::Zero(n);
    omega = 0;
    beta = 0;
    rk = b - A*xk;
    r0c = rk;
    pk = rk;
}

template<> double BiCGSTAB<std::complex<float>>::norm_sq(VectorXT vec){
    return (vec.dot(vec)).real();
}
template<> double BiCGSTAB<std::complex<double>>::norm_sq(VectorXT vec){
    return (vec.dot(vec)).real();
}
template<> double BiCGSTAB<double>::norm_sq(VectorXT vec){
    return vec.dot(vec);
}
template<> double BiCGSTAB<float>::norm_sq(VectorXT vec){
    return vec.dot(vec);
}
template<> double BiCGSTAB<int>::norm_sq(VectorXT vec){
    return vec.dot(vec);
}

template<typename T>
void BiCGSTAB<T>::solve(){
    int k = 0;
    int matvec_count = 0;
    std::ofstream logs("../output/BiCGSTAB_logs.csv", std::ios::out | std::ios::trunc);
    logs << "k,time,res_2norm_rel,matvec_count\n";
    auto start = std::chrono::high_resolution_clock::now();
    auto check = std::chrono::high_resolution_clock::now();
    while (k < n){
        vk = A*pk;
        matvec_count++;
        alpha = r0c.dot(rk)/r0c.dot(vk);
        sk = rk - alpha*vk;
        double sk_norm_sq_rel = norm_sq(sk)/norm_sq(r0c);      
        std::cout << "step: " << float(k) + 0.5 << ", (sk,sk)^1/2 / (r0,r0)^1/2 = " << std::sqrt(sk_norm_sq_rel) << "\n\n";
        check = std::chrono::high_resolution_clock::now();
        logs << float(k) + 0.5 << ',' 
             << std::chrono::duration_cast<std::chrono::microseconds>(check-start).count() << "," 
             << std::sqrt(sk_norm_sq_rel) << ',' 
             << matvec_count << '\n';
        if (sk_norm_sq_rel < epsilon*epsilon){
            xk += alpha * pk;
            break;
        }
        tk = A*sk;
        matvec_count++;
        omega = tk.dot(sk)/tk.dot(tk);
        xk += alpha * pk + omega * sk;
        rkp1 = sk - omega * tk;
        double rkp1_norm_sq_rel = norm_sq(rkp1)/norm_sq(r0c);
        std::cout << "step: " << k + 1 << ", (rk,rk)^1/2 / (r0c,r0c)^1/2 = " << std::sqrt(rkp1_norm_sq_rel) << "\n\n";// << "xk = \n" << xk << "\n\n\n\n";
        check = std::chrono::high_resolution_clock::now();
        logs << k + 1 << ',' 
             << std::chrono::duration_cast<std::chrono::microseconds>(check-start).count() << "," 
             << std::sqrt(rkp1_norm_sq_rel) << ',' 
             << matvec_count <<'\n';
        if (rkp1_norm_sq_rel < epsilon*epsilon)
            break;
        beta = r0c.dot(rkp1)/r0c.dot(rk)*alpha/omega;
        rk = rkp1;
        pk = rkp1 + beta * (pk - omega * vk);
        k+=1;
    }
    logs.close();
    write_binary("../output/BiCGSTAB_solution.dat", xk);
}