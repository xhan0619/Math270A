#include <cmath>
#include "Tools.h"
#include "HW1.h"
#include "runPolarBenchmark.h"
#include "iostream"

int main()
{	
// 	using namespace JIXIE;
// 	Eigen::Matrix<float, 2, 2> F;
// 	Eigen::Matrix<float, 2, 2> V;
// 	Eigen::Matrix<float, 2, 2> U;
// 	Eigen::Matrix<float, 2, 1> Sigma;
// 	float sigma1, sigma2;
// 	F << 2.000026941, 1.999972224,
// 1.999996662, 2.00002408
// ;

// 	bool det_V_negative = false;
// 	bool det_U_negative = false;

// 	// Make symmetric C = F^T * F
// 	Eigen::Matrix<float, 2, 2> C;
// 	C.noalias() = F.transpose() * F;

// 	// Jacobi Algorithm to get V and positive sigma1 and sigma2
// 	Jacobi(C, V, sigma1, sigma2);
// 	std::cout << "sigma2 is " << sigma2 << std::endl;
// 	sigma1 = sqrt(sigma1);
// 	sigma2 = sqrt(sigma2);

// 	Sort(V, sigma1, sigma2, det_V_negative); // Sort eigenvalues and V accordingly
// 	std::cout << "V is " << V << std::endl;
// 	// set up A for QR factorization
// 	Eigen::Matrix<float, 2, 2> A; 
// 	A.noalias() = F * V; 
// 	GivensRotation<float> r(A(0, 0), A(1, 0), 0, 1);

// 	// A becomes the R matrix after Givens
// 	r.rowRotation(A);

// 	// let U be the rotation r^T
// 	//r.transposeInPlace();
// 	r.fill(U);

// 	// If R (which is A now) has 2,2 entry negative,
// 	// U is [-s -c]
// 	if (A(1, 1) < 0){
// 		// flip the sign of second column of U
// 		U.col(1) = -U.col(1);
// 		det_U_negative = true;
// 	}

// 	std::cout << "U is " << U << std::endl;
// 	std::cout << "V is " << V << std::endl;
// 	std::cout << "Sigma is " << sigma1 << sigma2 << std::endl;
// 	enforceSignConvention(F, U, V, sigma1, sigma2, det_U_negative, det_V_negative);

// 	std::cout << "U is " << U << std::endl;
// 	std::cout << "V is " << V << std::endl;
// 	std::cout << "Sigma is " << sigma1 << sigma2 << std::endl;
// 	Sigma(0,0) = sigma1;
// 	Sigma(1,0) = sigma2;
// 	std::cout << U * Sigma.asDiagonal() * V.transpose() << std::endl;

	// std::cout << sigma1 << std::endl << sigma2 << std::endl;

	//  singularValueDecomposition(F, U, Sigma, V);
	//  std::cout << "U is " << U << std::endl;
	//  std::cout << "V is " << V << std::endl;
	//  std::cout << "Sigma is " << Sigma << std::endl;
	//  float error = (U * Sigma.asDiagonal() * V.transpose() - F).array().abs().maxCoeff();
	//  std::cout << error << std::endl;

	// std::cout << U * Sigma.asDiagonal() * V.transpose() << std::endl;



     bool run_benchmark = true;
     if (run_benchmark) runBenchmark();
}
