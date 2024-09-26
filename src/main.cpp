# include <iostream>
# include <Eigen/Dense>
# include <BiCGSTAB.h>
# include <BlBiCGSTAB.h>
# include <fbinio.h>

int main (int argc, char** argv){
// ------------Linear system AX = B-----------------      
    // Eigen::MatrixXcd A(3,3);
    // A << 2, 3, 5,
    //      3, 7, 4,
    //      1, 2, 2;  
    // std::cout << "A = \n" << A << "\n\n";
// -------------------------------------------------
// ------------------Many RHS Test------------------
    // Eigen::MatrixXcd b(3, 2);
    // b << 10, 0,
    //       3, 1,
    //       3, 0;
    // std::cout << "b = \n" << b << "\n\n";
    // Eigen::MatrixXcd x0(3, 2);
    // x0 << 0, 0,
    //       0, 0,
    //       0, 0;
    // std::cout << "x0 = \n" << x0 << "\n\n";
// ------------------------------------------------
//-------------------One RHS Test------------------
    // Eigen::VectorXd b(3);
    // b << 0, 1, 0;
    // std::cout << b << "\n\n";
    // Eigen::VectorXd x0(3);
    // x0 << 0, 0 ,0; 
    // std::cout << x0 << "\n\n";

    // BiCGSTAB bcg(A, b, x0, 0.001);
    // bcg.solve();
    // std::cout << "Solution = \n" << bcg.xk << "\n\n";

// ------------------------------------------------
//--------------BlBiCGSTAB Solve-------------------
    // BlBiCGSTAB<Eigen::MatrixXcd, Eigen::VectorXcd> bcg(A, b, x0, 0.001);
    // bcg.solve();
    // std::cout << "Solution = \n" << bcg.Xk << "\n\n";
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
//-----------------------fbinio test-----------------
// Eigen::MatrixXd mat(10, 10);
// Eigen::MatrixXd matwr = Eigen::MatrixXd::Random(10, 10);
// write_binary<Eigen::MatrixXd>("matwr.dat", matwr);
// read_binary<Eigen::MatrixXd>("matwr.dat", mat);
// std::cout << mat << "\n\n";
//------------------------------------------------

//-----------------BiCGSTAB combat task----------------------
Eigen::MatrixXcf rhs(3, 3); // "3" is arbitrary, in read_binary() function it will
Eigen::MatrixXcf A(3,3);    // be resized

read_binary("/home/starman/rhs_alm_722.dat", rhs);
read_binary("/home/starman/mat_alm_full.dat", A);
Eigen::VectorXcf b = rhs.col(0);
std::cout << "A is " << A.rows() << "x" << A.cols() << "\n\n";
std::cout << "b is " << b.rows() << "x" << b.cols() << "\n\n";
Eigen::VectorXcf x0 = Eigen::VectorXcf::Zero(b.rows()); 
BiCGSTAB bcg(A, b, x0, 0.001);
bcg.solve();
//-----------------BlBiCG combat task-------------------------
// Eigen::MatrixXcf rhs(3, 3); // "3" is arbitrary, in read_binary() function it will
// Eigen::MatrixXcf A(3,3);    // be resized

// read_binary("/home/starman/rhs_alm_722.dat", rhs);
// read_binary("/home/starman/mat_alm_full.dat", A);

// Eigen::MatrixXcf B = rhs.leftCols(1);

// std::cout << "A is " << A.rows() << "x" << A.cols() << "\n\n";
// std::cout << "B is " << B.rows() << "x" << B.cols() << "\n\n";

// Eigen::MatrixXcf X0 = Eigen::MatrixXcf::Zero(B.rows(), B.cols());

// BlBiCGSTAB<Eigen::MatrixXcf, Eigen::VectorXcf> bbcg(A, B, X0, 0.001);
// bbcg.solve();
//------------------------------------------------------------

return 0;
}