#ifndef JIXIE_HW1_H
#define JIXIE_HW1_H

#include "Tools.h"
#include <iostream>
#include <cmath>
#include <vector>
/**
SVD based on implicit QR with Wilkinson Shift
*/

namespace JIXIE {

/**
    Class for givens rotation (copied from ImplicitQRSVD.h).
    Row rotation G*A corresponds to something like
    c -s  0
    ( s  c  0 ) A
    0  0  1
    Column rotation A G' corresponds to something like
    c -s  0
    A ( s  c  0 )
    0  0  1

    c and s are always computed so that
    ( c -s ) ( a )  =  ( * )
    s  c     b       ( 0 )

    Assume rowi<rowk.
    */
template <class T>
class GivensRotation {
public:
    int rowi;
    int rowk;
    T c;
    T s;

    inline GivensRotation(int rowi_in, int rowk_in)
        : rowi(rowi_in)
        , rowk(rowk_in)
        , c(1)
        , s(0)
    {
    }

    inline GivensRotation(T a, T b, int rowi_in, int rowk_in)
        : rowi(rowi_in)
        , rowk(rowk_in)
    {
        compute(a, b);
    }


    ~GivensRotation() {}

    inline void transposeInPlace()
    {
        s = -s;
    }

    /**
        Compute c and s from a and b so that
        ( c -s ) ( a )  =  ( * )
        s  c     b       ( 0 )
        */
    inline void compute(const T a, const T b)
    {
        using std::sqrt;

        T d = a * a + b * b;
        c = 1;
        s = 0;
        if (d != 0) {
            // T t = 1 / sqrt(d);
            T t = JIXIE::MATH_TOOLS::rsqrt(d);
            c = a * t;
            s = -b * t;
        }
    }

    /**
        This function computes c and s so that
        ( c -s ) ( a )  =  ( 0 )
        s  c     b       ( * )
        */
    inline void computeUnconventional(const T a, const T b)
    {
        using std::sqrt;

        T d = a * a + b * b;
        c = 0;
        s = 1;
        if (d != 0) {
            // T t = 1 / sqrt(d);
            T t = JIXIE::MATH_TOOLS::rsqrt(d);
            s = a * t;
            c = b * t;
        }
    }

    /**
      Fill the R with the entries of this rotation
        */
    template <class MatrixType>
    inline void fill(const MatrixType& R) const
    {
        MatrixType& A = const_cast<MatrixType&>(R);
        A = MatrixType::Identity();
        A(rowi, rowi) = c;
        A(rowk, rowi) = -s;
        A(rowi, rowk) = s;
        A(rowk, rowk) = c;
    }

    /**
        This function does something like
        c -s  0
        ( s  c  0 ) A -> A
        0  0  1
        It only affects row i and row k of A.
        */
    template <class MatrixType>
    inline void rowRotation(MatrixType& A) const
    {
        for (int j = 0; j < MatrixType::ColsAtCompileTime; j++) {
            T tau1 = A(rowi, j);
            T tau2 = A(rowk, j);
            A(rowi, j) = c * tau1 - s * tau2;
            A(rowk, j) = s * tau1 + c * tau2;
        }
    }

    /**
        This function does something like
        c  s  0
        A ( -s  c  0 )  -> A
        0  0  1
        It only affects column i and column k of A.
        */
    template <class MatrixType>
    inline void columnRotation(MatrixType& A) const
    {
        for (int j = 0; j < MatrixType::RowsAtCompileTime; j++) {
            T tau1 = A(j, rowi);
            T tau2 = A(j, rowk);
            A(j, rowi) = c * tau1 - s * tau2;
            A(j, rowk) = s * tau1 + c * tau2;
        }
    }

    /**
      Multiply givens must be for same row and column
      **/
    inline void operator*=(const GivensRotation<T>& A)
    {
        T new_c = c * A.c - s * A.s;
        T new_s = s * A.c + c * A.s;
        c = new_c;
        s = new_s;
    }

    /**
      Multiply givens must be for same row and column
      **/
    inline GivensRotation<T> operator*(const GivensRotation<T>& A) const
    {
        GivensRotation<T> r(*this);
        r *= A;
        return r;
    }
};

/**
       \brief 2X2 Jacobi routine.
       \param[in] C matrix positive semidefinite.
       \param[out] V Robustly a rotation matrix.
       \param[out] sigma1 positive eigenvalue with larger magnitude.
       \param[out] sigma2 eigenvalue with smaller magnitude.

       Whole matrix V is stored
       Guarantees negative sign is on the small magnitude singular value.
*/
template <class T>
inline void Jacobi(const Eigen::Matrix<T, 2, 2>& C,
    Eigen::Matrix<T, 2, 2>& V,
    T& sigma1,
    T& sigma2)
{
	using std::sqrt;
	// if already diagonal, take V = I_2
	if (std::abs(C(0, 1)) == 0){
		V = Eigen::Matrix<T, 2, 2>::Identity();
		sigma1 = C(0, 0);
		sigma2 = C(1, 1);
		return;
	}

	T tau =  (C(1, 1) - C(0, 0)) / (2 * C(0, 1));
	T t2;
	Eigen::Matrix<T, 2, 2> Sigma; // Eigenvalues

	// Assign t2 to avoid numerical inprecision
	// TODO: check the conditions
	if (tau > 0){
		//if(true){
			t2 = (-1) / (tau + sqrt(1 + tau * tau));
	//	}
	//	else {t2 = C(1, 0) / (sqrt(C(1, 0) * C(1, 0) + (C(1, 1) - C(0, 0)) * (C(1, 1) - C(0, 0)) / 4) + (C(1, 1) - C(0, 0)) / 2); }
	}
	else {
//		if (true){
			t2 = 1 / (sqrt(1 + tau * tau) - tau);
//		}
//		else {t2 = C(1, 0) / (sqrt(C(1, 0) * C(1, 0) + (C(1, 1) - C(0, 0)) * (C(1, 1) - C(0, 0)) / 4) - (C(1, 1) - C(0, 0)) / 2); }
	}

	T c2 = JIXIE::MATH_TOOLS::rsqrt(1 + t2 * t2); // define c2
	T s2 = t2 * c2; // define s2
	V << c2, (-1)*s2, 
		 s2, c2; // update V

	Sigma.noalias() = V.transpose() * C * V; // update C  
	sigma1 = std::max(T(Sigma(0, 0)), T(0));
	sigma2 = std::max(T(Sigma(1, 1)), T(0));


    // sigma1 = std::max(T(pow(c2, 2) * C(0, 0) + 2. * c2 * s2 * C(1, 0) + pow(s2, 2) * C(1, 1)), T(0));
    // sigma2 = std::max(T(pow(s2, 2) * C(0, 0) - 2. * c2 * s2 * C(1, 0) + pow(c2, 2) * C(1, 1)), T(0));

};

/**
Helper function to sort the singular values in decreasing magnitude
and permute columns accordingly
*/
template <class T>
inline void Sort(Eigen::Matrix<T, 2, 2>& V,
    T& sigma1,
    T& sigma2,
    bool& det_V_negative)
{
	if (sigma2 > sigma1){
		// swap sigma1 and sigma2
		T temp = sigma2;
		sigma2 = sigma1;
		sigma1 = temp;

		// swap the columns of V
		V.col(0).swap(V.col(1));

		// mark the flag as TRUE
		det_V_negative = true;
	}
}

template <class T>
inline void enforceSignConvention(const Eigen::Matrix<T, 2, 2>& F,
	Eigen::Matrix<T, 2, 2>& U,
	Eigen::Matrix<T, 2, 2>& V,
    T& sigma1,
    T& sigma2,
    bool& det_U_negative,
    bool& det_V_negative)
{
	T det = F(0, 0) * F(1, 1) - F(0, 1) * F(1, 0);
	if (det < 0){
		sigma2 = -sigma2;
		if (det_U_negative){
			// flip the sign of second column of U
			U.col(1) = -U.col(1);
		}
		else {
			// flip the sign of second column of V
			V.col(1) = -V.col(1);
		}
	}
	else if (det > 0){
		if (det_U_negative && det_V_negative){
			// flip the sign of second column of U
			U.col(1) = -U.col(1);
			// flip the sign of second column of V
			V.col(1) = -V.col(1);
		}
	}
	else {
		if (det_V_negative){
			// flip the sign of second column of V
			V.col(0) = -V.col(0);

			// if only one of U, V flipped, flip sigma1
			if (!det_U_negative){
				sigma1 = -sigma1;
			}
		}
		if (det_U_negative){
			// flip the sign of second column of U
			U.col(0) = -U.col(0);
			// if only one of U, V flipped, flip sigma1
			if (!det_V_negative){
				sigma1 = -sigma1;
			}
		}
	}
}

/**
   \brief 2x2 SVD (singular value decomposition) A=USV'
   \param[in] A Input matrix.
   \param[out] U Robustly a rotation matrix.
   \param[out] Sigma Vector of singular values sorted with decreasing magnitude. The second one can be negative.
   \param[out] V Robustly a rotation matrix.
*/
template <class T>
inline int singularValueDecomposition(
	const Eigen::Matrix<T, 2, 2>& F,
    Eigen::Matrix<T, 2, 2>& U,
    Eigen::Matrix<T, 2, 1>& Sigma,
    Eigen::Matrix<T, 2, 2>& V)
{
	using std::sqrt;
	// initialize the flag for det of V and U
	bool det_V_negative = false;
	bool det_U_negative = false;

	// Make symmetric C = F^T * F
	Eigen::Matrix<T, 2, 2> C;
	C.noalias() = F.transpose() * F;

	T sigma1, sigma2 = 0;
	// Jacobi Algorithm to get V and positive sigma1 and sigma2
	Jacobi(C, V, sigma1, sigma2);

	sigma1 = sqrt(sigma1);
	sigma2 = sqrt(sigma2);

	Sort(V, sigma1, sigma2, det_V_negative); // Sort eigenvalues and V accordingly

	// set up A for QR factorization
	Eigen::Matrix<T, 2, 2> A; 
	A.noalias() = F * V; 
	GivensRotation<T> r(A(0, 0), A(1, 0), 0, 1);

	// A becomes the R matrix after Givens
	r.rowRotation(A);

	// let U be the rotation r^T
	//r.transposeInPlace();
	r.fill(U);

	// If R (which is A now) has 2,2 entry negative,
	// U is [-s -c]
	if (A(1, 1) < 0){
		// flip the sign of second column of U
		U.col(1) = -U.col(1);
		det_U_negative = true;
	}

	// At this point we shoule have 
	// F = tilda(U) * tilda(Sigma) * tilda(V)^T
	// with tilda(Sigma) positive

	enforceSignConvention(F, U, V, sigma1, sigma2, det_U_negative, det_V_negative);
	Sigma(0,0) = sigma1;
	Sigma(1,0) = sigma2;
	//Sigma << sigma1, sigma2;
	return 0;
}

template <class T>
inline int polarDecomposition(
	const Eigen::Matrix<T, 3, 3>& F,
    Eigen::Matrix<T, 3, 3>& R,
    Eigen::Matrix<T, 3, 3>& S,
    T tol = 128 * std::numeric_limits<T>::epsilon())
{
	R = Eigen::Matrix<T, 3, 3>::Identity();
	S = F;
	int it = 0;
	int max_it = 500;

    GivensRotation<T> r01(0, 1);
    GivensRotation<T> r02(0, 2);
    GivensRotation<T> r12(1, 2);
    GivensRotation<T> vec[] = {r01, r02, r12};

    //while ((it < max_it) && (std::max({std::abs(S(0, 1) - S(1, 0)), std::abs(S(0, 2) - S(2, 0)), std::abs(S(1, 2) - S(2, 1))}) > tol)){
	while ((it < max_it) && ((std::abs(S(0, 1) - S(1, 0)) > tol) || (std::abs(S(0, 2) - S(2, 0)) > tol) || (std::abs(S(1, 2) - S(2, 1)) > tol))){
		for (int j = 1; j < 3; j++){
			for (int i = 0; i < j; i++){
                vec[i + j - 1].compute(S(i, i) + S(j, j), S(j, i) - S(i, j));
                vec[i + j - 1].columnRotation(R);
                vec[i + j - 1].rowRotation(S);

				// GivensRotation<T> r(S(i, i) + S(j, j), S(j, i) - S(i, j), i, j);
				// r.columnRotation(R);
				// r.rowRotation(S);
			}
		}
		++it;
	}
	return it;
}
}
#endif