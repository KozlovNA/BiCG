# include <iostream>
# include <Eigen/Dense>
# include <BiCGSTAB.h>
# include <BlBiCGSTAB.h>

int main (int argc, char** argv){
// ------------Linear system AX = B-----------------      
    Eigen::MatrixXd A(3,3);
    A << 2, 3, 5,
         3, 7, 4,
         1, 2, 2;  
    std::cout << "A = \n" << A << "\n\n";
// -------------------------------------------------
// ------------------Many RHS Test------------------
    Eigen::MatrixXd b(3, 2);
    b << 10, 0,
          3, 1,
          3, 0;
    std::cout << "b = \n" << b << "\n\n";
    Eigen::MatrixXd x0(3, 2);
    x0 << 0, 0,
          0, 0,
          0, 0;
    std::cout << "x0 = \n" << x0 << "\n\n";
// ------------------------------------------------
//-------------------One RHS Test------------------
//      BlBiCGSTAB bcg(A, b, x0, 0.001);          
//     Eigen::Vector3d b(0, 1, 0);
//     std::cout << b << "\n\n";
//     Eigen::Vector3d x0(0, 0, 0);
// ------------------------------------------------
//--------------BlBiCGSTAB Solve-------------------
    BlBiCGSTAB bcg(A, b, x0, 0.001);
    bcg.solve();
    std::cout << "Solution = \n" << bcg.Xk << "\n\n";
//-------------------------------------------------    
// ----------------L_solve Test--------------------
//     Eigen::VectorXd b(3);
//     b << 1,
//          2,
//          6;
//     std::cout << b << "\n\n";
//     Eigen::MatrixXd L(3, 3);
//     L << 1, 0, 0,
//          1, 1, 0,
//          2, 3, 1;
//     std::cout << L << "\n\n";
//     Eigen::VectorXd x(3);     
//     x = L_solve(L, b);
//     std::cout << x << "\n\n";
// ---------------------------------------------- 
// ----------------U_solve Test-----------------
      // Eigen::VectorXd b(3);
      // b <<  6,
      //       2,
      //       1;
      // std::cout << b << "\n\n";
      // Eigen::MatrixXd U(3, 3);
      // U << 1, 2, 3,
      // 0, 1, 1,
      // 0, 0, 1;
      // std::cout << U << "\n\n";
      // Eigen::VectorXd x(3);     
      // x = U_solve(U, b);
      // std::cout << x << "\n\n";
// ----------------------------------------------
// ------------------LU_solve Test----------------
//     Eigen::Vector3d tmp(0, 1, 0);
//     Eigen::MatrixXd rhs = tmp;
//     std::cout << rhs << "\n\n";
//     Eigen::MatrixXd sol = LU_solve(A, rhs);
//     std::cout << sol << "\n\n";     
return 0;
}