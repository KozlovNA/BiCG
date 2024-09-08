# include <iostream>
# include <Eigen/Dense>
# include <BiCGSTAB.h>
# include <BlBiCGSTAB.h>

int main (int argc, char** argv){
    Eigen::MatrixXd A(3,3);
    A << 2, 3, 5,
         3, 7, 4,
         1, 2, 2;  
    std::cout << A << "\n\n";
    Eigen::Vector3d b(10, 3, 3);
    std::cout << b << "\n\n";
    Eigen::Vector3d x0(0, 0, 0);
    std::cout << x0 << "\n\n";
    BlBiCGSTAB bcg(A, b, x0, 2);
    bcg.solve(); 
    std::cout << bcg.Xk << "\n\n";
    return 0;
}