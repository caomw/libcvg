//
// sqpnp.cpp
//
// George Terzakis (terzakig-at-hotmail-dot-com), September 2020
// 
// Implementation of SQPnP as described in the paper:
//
// "A Consistently Fast and Globally Optimal Solution to the Perspective-n-Point Problem" by G. Terzakis and M. Lourakis
//  	 a) Paper: 	   http://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf 
//       b) Supplementary: https://www.ecva.net/papers/eccv_2020/papers_ECCV/papers/123460460.pdf


#pragma warning(disable:4996)

#include <cmath>
#include <iostream>
#include <fstream>

#include <boost/algorithm/string.hpp>

#include <Eigen/Eigen>

#include "sqpnp.h"


namespace cvg
{

	const double SolverParameters::DEFAULT_RANK_TOLERANCE = 1e-7;
	const double SolverParameters::DEFAULT_SQP_SQUARED_TOLERANCE = 1e-10;
	const double SolverParameters::DEFAULT_SQP_DET_THRESHOLD = 1.001;
	const NearestRotationMethod SolverParameters::DEFAULT_NEAREST_ROTATION_METHOD = NearestRotationMethod::FOAM;
	const double SolverParameters::DEFAULT_ORTHOGONALITY_SQUARED_ERROR_THRESHOLD = 1e-8;
	const double SolverParameters::DEFAULT_EQUAL_VECTORS_SQUARED_DIFF = 1e-10;
	const double SolverParameters::DEFAULT_EQUAL_SQUARED_ERRORS_DIFF = 1e-6;
	const double SolverParameters::DEFAULT_POINT_VARIANCE_THRESHOLD = 1e-5;


	const double SqPnPSolver::SQRT3 = std::sqrt(3);


	void SqPnPSolver::HandleSolution(SQPSolution& solution, double& min_sq_error)
	{
		if (TestPositiveDepth(solution))
		{

			solution.sq_error = (Omega_ * solution.r_hat).dot(solution.r_hat);
			if (std::fabs(min_sq_error - solution.sq_error) > parameters_.equal_squared_errors_diff)
			{
				if (min_sq_error > solution.sq_error)
				{
					min_sq_error = solution.sq_error;
					solutions_[0] = solution;
					num_solutions_ = 1;
				}
			}
			else // look for a solution that's almost equal to this
			{
				bool found = false;
				for (int i = 0; i < num_solutions_; i++)
				{
					if ((solutions_[i].r_hat - solution.r_hat).squaredNorm() < parameters_.equal_vectors_squared_diff)
					{
						if (solutions_[i].sq_error > solution.sq_error)
						{
							solutions_[i] = solution;
						}
						found = true;
						break;
					}
				}

				if (!found)
				{
					solutions_[num_solutions_++] = solution;
				}

				if (min_sq_error > solution.sq_error)
				{
					min_sq_error = solution.sq_error;
				}
			}
		}
	}


	// Constructor (initializes Omega and P and U, s, i.e. the decomposition of Omega)
	template <class Point3D, class Projection2D>
	void SqPnPSolver::indirectSolver(
		const std::vector<Point3D>& _3dpoints,
		const std::vector<Projection2D>& _projections,
		const SolverParameters& _parameters /*= SolverParameters()*/)
	{


		this->parameters_ = _parameters;

		if (_3dpoints.size() != _projections.size() ||
			_3dpoints.size() < 3 || _projections.size() < 3)
		{
			flag_valid_ = false;
			return;
		}

		flag_valid_ = true;
		num_null_vectors_ = -1; // set to -1 in case we never make it to the decomposition of Omega
		Omega_ = Eigen::Matrix<double, 9, 9>::Zero();
		const size_t n = _3dpoints.size();
		double sum_x = 0,
			sum_y = 0,
			sum_x2_plus_y2 = 0;
		double sum_X = 0,
			sum_Y = 0,
			sum_Z = 0;

		Eigen::Matrix<double, 3, 9> QA = Eigen::Matrix<double, 3, 9>::Zero();  // Sum( Qi*Ai )

		for (size_t i = 0; i < n; i++)
		{
			points_.push_back(_3dpoints.at(i));
			projections_.push_back(_projections.at(i));

			double x = projections_.rbegin()->vector[0],
				y = projections_.rbegin()->vector[1],
				sq_norm_m = projections_.rbegin()->vector.squaredNorm();
			sum_x += x;
			sum_y += y;
			sum_x2_plus_y2 += sq_norm_m;

			double X = points_.rbegin()->vector[0],
				Y = points_.rbegin()->vector[1],
				Z = points_.rbegin()->vector[2];
			sum_X += X;
			sum_Y += Y;
			sum_Z += Z;

			// Accumulate Omega by kronecker( Qi, Mi*Mi' ) = Q + A'*Qi*Ai. NOTE: Skipping block (3:5, 3:5) because its same as (0:2, 0:2)
			double X2 = X * X, XY = X * Y, XZ = X * Z, Y2 = Y * Y, YZ = Y * Z, Z2 = Z * Z;
			// a. Block (0:2, 0:2) populated by Mi*Mi'. NOTE: Only upper triangle
			Omega_(0, 0) += X2;
			Omega_(0, 1) += XY;
			Omega_(0, 2) += XZ;
			Omega_(1, 1) += Y2;
			Omega_(1, 2) += YZ;
			Omega_(2, 2) += Z2;

			// b. Block (0:2, 6:8) populated by -x*Mi*Mi'. NOTE: Only upper triangle
			Omega_(0, 6) += -x * X2; Omega_(0, 7) += -x * XY; Omega_(0, 8) += -x * XZ;
			Omega_(1, 7) += -x * Y2; Omega_(1, 8) += -x * YZ;
			Omega_(2, 8) += -x * Z2;
			// c. Block (3:5, 6:8) populated by -y*Mi*Mi'. NOTE: Only upper triangle
			Omega_(3, 6) += -y * X2; Omega_(3, 7) += -y * XY; Omega_(3, 8) += -y * XZ;
			Omega_(4, 7) += -y * Y2; Omega_(4, 8) += -y * YZ;
			Omega_(5, 8) += -y * Z2;
			// d. Block (6:8, 6:8) populated by (x^2+y^2)*Mi*Mi'. NOTE: Only upper triangle
			Omega_(6, 6) += sq_norm_m * X2; Omega_(6, 7) += sq_norm_m * XY; Omega_(6, 8) += sq_norm_m * XZ;
			Omega_(7, 7) += sq_norm_m * Y2; Omega_(7, 8) += sq_norm_m * YZ;
			Omega_(8, 8) += sq_norm_m * Z2;

			// Accumulating Qi*Ai in QA
			QA(0, 0) += X; QA(0, 1) += Y; QA(0, 2) += Z; 	QA(0, 6) += -x * X; QA(0, 7) += -x * Y; QA(0, 8) += -x * Z;
			QA(1, 3) += X; QA(1, 4) += Y; QA(1, 5) += Z; 	QA(1, 6) += -y * X; QA(1, 7) += -y * Y; QA(1, 8) += -y * Z;

			QA(2, 0) += -x * X; QA(2, 1) += -x * Y; QA(2, 2) += -x * Z; 	QA(2, 3) += -y * X; QA(2, 4) += -y * Y; QA(2, 5) += -y * Z;
			QA(2, 6) += sq_norm_m * X; QA(2, 7) += sq_norm_m * Y; QA(2, 8) += sq_norm_m * Z;

		}

		// Fill-in lower triangles of off-diagonal blocks (0:2, 6:8), (3:5, 6:8) and (6:8, 6:8)
		Omega_(1, 6) = Omega_(0, 7); Omega_(2, 6) = Omega_(0, 8); Omega_(2, 7) = Omega_(1, 8);
		Omega_(4, 6) = Omega_(3, 7); Omega_(5, 6) = Omega_(3, 8); Omega_(5, 7) = Omega_(4, 8);
		Omega_(7, 6) = Omega_(6, 7); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);

		// Fill-in upper triangle of block (3:5, 3:5)
		Omega_(3, 3) = Omega_(0, 0); Omega_(3, 4) = Omega_(0, 1); Omega_(3, 5) = Omega_(0, 2);
		Omega_(4, 4) = Omega_(1, 1); Omega_(4, 5) = Omega_(1, 2);
		Omega_(5, 5) = Omega_(2, 2);
		// Fill lower triangle of Omega
		Omega_(1, 0) = Omega_(0, 1);
		Omega_(2, 0) = Omega_(0, 2); Omega_(2, 1) = Omega_(1, 2);
		Omega_(3, 0) = Omega_(0, 3); Omega_(3, 1) = Omega_(1, 3); Omega_(3, 2) = Omega_(2, 3);
		Omega_(4, 0) = Omega_(0, 4); Omega_(4, 1) = Omega_(1, 4); Omega_(4, 2) = Omega_(2, 4); Omega_(4, 3) = Omega_(3, 4);
		Omega_(5, 0) = Omega_(0, 5); Omega_(5, 1) = Omega_(1, 5); Omega_(5, 2) = Omega_(2, 5); Omega_(5, 3) = Omega_(3, 5); Omega_(5, 4) = Omega_(4, 5);
		Omega_(6, 0) = Omega_(0, 6); Omega_(6, 1) = Omega_(1, 6); Omega_(6, 2) = Omega_(2, 6); Omega_(6, 3) = Omega_(3, 6); Omega_(6, 4) = Omega_(4, 6); Omega_(6, 5) = Omega_(5, 6);
		Omega_(7, 0) = Omega_(0, 7); Omega_(7, 1) = Omega_(1, 7); Omega_(7, 2) = Omega_(2, 7); Omega_(7, 3) = Omega_(3, 7); Omega_(7, 4) = Omega_(4, 7); Omega_(7, 5) = Omega_(5, 7); Omega_(7, 6) = Omega_(6, 7);
		Omega_(8, 0) = Omega_(0, 8); Omega_(8, 1) = Omega_(1, 8); Omega_(8, 2) = Omega_(2, 8); Omega_(8, 3) = Omega_(3, 8); Omega_(8, 4) = Omega_(4, 8); Omega_(8, 5) = Omega_(5, 8); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);

		// Q = Sum( Qi ) = Sum( [ 1, 0, -x; 0, 1, -y; -x, -y, x^2 + y^2] )
		Eigen::Matrix<double, 3, 3> Q;
		Q(0, 0) = n; Q(0, 1) = 0; Q(0, 2) = -sum_x;
		Q(1, 0) = 0; Q(1, 1) = n; Q(1, 2) = -sum_y;
		Q(2, 0) = -sum_x; Q(2, 1) = -sum_y; Q(2, 2) = sum_x2_plus_y2;

		// Qinv = inv( Q ) = inv( Sum( Qi) )
		double inv_n = 1.0 / n;
		double detQ = n * (n*sum_x2_plus_y2 - sum_y * sum_y - sum_x * sum_x);
		double point_coordinate_variance = detQ * inv_n * inv_n * inv_n;
		if (point_coordinate_variance < parameters_.point_variance_threshold)
		{
			flag_valid_ = false;
			return;
		}

		Eigen::Matrix<double, 3, 3> Qinv;
		InvertSymmetric3x3(Q, Qinv);

		// Compute P = -inv( Sum(Qi) ) * Sum( Qi*Ai ) = -Qinv * QA
		P_ = -Qinv * QA;
		// Complete Omega (i.e., Omega = Sum(A'*Qi*A') + Sum(Qi*Ai)'*P = Sum(A'*Qi*A') + Sum(Qi*Ai)'*inv(Sum(Qi))*Sum( Qi*Ai) 
		Omega_ += QA.transpose()*P_;

		// Finally, decompose Omega
		Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(Omega_, Eigen::ComputeFullU);

		U_ = svd.matrixU();
		s_ = svd.singularValues();

		// Find dimension of null space
		while (s_[7 - num_null_vectors_] < _parameters.rank_tolerance) num_null_vectors_++;
		// Dimension of null space of Omega must be <= 6
		if (++num_null_vectors_ > 6)
		{
			flag_valid_ = false;
		}

		// Point mean 
		point_mean_ << sum_X * inv_n, sum_Y*inv_n, sum_Z*inv_n;

		// Assign nearest rotation method
		if (parameters_.nearest_rotation_method == NearestRotationMethod::FOAM)
		{
			NearestRotationMatrix = NearestRotationMatrix_FOAM;
		}
		else // if ( parameters_.nearest_rotation_method == NearestRotationMethod::SVD )
		{
			NearestRotationMatrix = NearestRotationMatrix_SVD;
		}
	}


	void SqPnPSolver::ComputeRT(
		std::vector<_Point> XXw,
		std::vector<_Projection> xxn,
		Eigen::Matrix3d &R,
		Eigen::Vector3d &t,
		const SolverParameters& _parameters /*= SolverParameters()*/)
	{

		//以下是构造函数中的内容
		std::vector<_Point>_3dpoints = XXw;
		std::vector<_Projection> _projections = xxn;

		if (_3dpoints.size() != _projections.size() ||
			_3dpoints.size() < 3 || _projections.size() < 3)
		{
			flag_valid_ = false;
			return;
		}

		flag_valid_ = true;
		num_null_vectors_ = -1; // set to -1 in case we never make it to the decomposition of Omega
		Omega_ = Eigen::Matrix<double, 9, 9>::Zero();
		const size_t n = _3dpoints.size();
		double sum_x = 0,
			sum_y = 0,
			sum_x2_plus_y2 = 0;

		double sum_X = 0,
			sum_Y = 0,
			sum_Z = 0;

		Eigen::Matrix<double, 3, 9> QA = Eigen::Matrix<double, 3, 9>::Zero();  // Sum( Qi*Ai )

		for (size_t i = 0; i < n; i++)
		{
			points_.push_back(_3dpoints.at(i));
			projections_.push_back(_projections.at(i));

			double x = projections_.rbegin()->vector[0],
				y = projections_.rbegin()->vector[1],
				sq_norm_m = projections_.rbegin()->vector.squaredNorm();
			sum_x += x;
			sum_y += y;
			sum_x2_plus_y2 += sq_norm_m;

			double X = points_.rbegin()->vector[0],
				Y = points_.rbegin()->vector[1],
				Z = points_.rbegin()->vector[2];
			sum_X += X;
			sum_Y += Y;
			sum_Z += Z;

			// Accumulate Omega by kronecker( Qi, Mi*Mi' ) = Q + A'*Qi*Ai. NOTE: Skipping block (3:5, 3:5) because its same as (0:2, 0:2)
			double X2 = X * X, XY = X * Y, XZ = X * Z, Y2 = Y * Y, YZ = Y * Z, Z2 = Z * Z;
			// a. Block (0:2, 0:2) populated by Mi*Mi'. NOTE: Only upper triangle
			Omega_(0, 0) += X2;
			Omega_(0, 1) += XY;
			Omega_(0, 2) += XZ;
			Omega_(1, 1) += Y2;
			Omega_(1, 2) += YZ;
			Omega_(2, 2) += Z2;

			// b. Block (0:2, 6:8) populated by -x*Mi*Mi'. NOTE: Only upper triangle
			Omega_(0, 6) += -x * X2; Omega_(0, 7) += -x * XY; Omega_(0, 8) += -x * XZ;
			Omega_(1, 7) += -x * Y2; Omega_(1, 8) += -x * YZ;
			Omega_(2, 8) += -x * Z2;
			// c. Block (3:5, 6:8) populated by -y*Mi*Mi'. NOTE: Only upper triangle
			Omega_(3, 6) += -y * X2; Omega_(3, 7) += -y * XY; Omega_(3, 8) += -y * XZ;
			Omega_(4, 7) += -y * Y2; Omega_(4, 8) += -y * YZ;
			Omega_(5, 8) += -y * Z2;
			// d. Block (6:8, 6:8) populated by (x^2+y^2)*Mi*Mi'. NOTE: Only upper triangle
			Omega_(6, 6) += sq_norm_m * X2; Omega_(6, 7) += sq_norm_m * XY; Omega_(6, 8) += sq_norm_m * XZ;
			Omega_(7, 7) += sq_norm_m * Y2; Omega_(7, 8) += sq_norm_m * YZ;
			Omega_(8, 8) += sq_norm_m * Z2;

			// Accumulating Qi*Ai in QA
			QA(0, 0) += X; QA(0, 1) += Y; QA(0, 2) += Z; 	QA(0, 6) += -x * X; QA(0, 7) += -x * Y; QA(0, 8) += -x * Z;
			QA(1, 3) += X; QA(1, 4) += Y; QA(1, 5) += Z; 	QA(1, 6) += -y * X; QA(1, 7) += -y * Y; QA(1, 8) += -y * Z;

			QA(2, 0) += -x * X; QA(2, 1) += -x * Y; QA(2, 2) += -x * Z; 	QA(2, 3) += -y * X; QA(2, 4) += -y * Y; QA(2, 5) += -y * Z;
			QA(2, 6) += sq_norm_m * X; QA(2, 7) += sq_norm_m * Y; QA(2, 8) += sq_norm_m * Z;

		}

		// Fill-in lower triangles of off-diagonal blocks (0:2, 6:8), (3:5, 6:8) and (6:8, 6:8)
		Omega_(1, 6) = Omega_(0, 7); Omega_(2, 6) = Omega_(0, 8); Omega_(2, 7) = Omega_(1, 8);
		Omega_(4, 6) = Omega_(3, 7); Omega_(5, 6) = Omega_(3, 8); Omega_(5, 7) = Omega_(4, 8);
		Omega_(7, 6) = Omega_(6, 7); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);

		// Fill-in upper triangle of block (3:5, 3:5)
		Omega_(3, 3) = Omega_(0, 0); Omega_(3, 4) = Omega_(0, 1); Omega_(3, 5) = Omega_(0, 2);
		Omega_(4, 4) = Omega_(1, 1); Omega_(4, 5) = Omega_(1, 2);
		Omega_(5, 5) = Omega_(2, 2);
		// Fill lower triangle of Omega
		Omega_(1, 0) = Omega_(0, 1);
		Omega_(2, 0) = Omega_(0, 2); Omega_(2, 1) = Omega_(1, 2);
		Omega_(3, 0) = Omega_(0, 3); Omega_(3, 1) = Omega_(1, 3); Omega_(3, 2) = Omega_(2, 3);
		Omega_(4, 0) = Omega_(0, 4); Omega_(4, 1) = Omega_(1, 4); Omega_(4, 2) = Omega_(2, 4); Omega_(4, 3) = Omega_(3, 4);
		Omega_(5, 0) = Omega_(0, 5); Omega_(5, 1) = Omega_(1, 5); Omega_(5, 2) = Omega_(2, 5); Omega_(5, 3) = Omega_(3, 5); Omega_(5, 4) = Omega_(4, 5);
		Omega_(6, 0) = Omega_(0, 6); Omega_(6, 1) = Omega_(1, 6); Omega_(6, 2) = Omega_(2, 6); Omega_(6, 3) = Omega_(3, 6); Omega_(6, 4) = Omega_(4, 6); Omega_(6, 5) = Omega_(5, 6);
		Omega_(7, 0) = Omega_(0, 7); Omega_(7, 1) = Omega_(1, 7); Omega_(7, 2) = Omega_(2, 7); Omega_(7, 3) = Omega_(3, 7); Omega_(7, 4) = Omega_(4, 7); Omega_(7, 5) = Omega_(5, 7); Omega_(7, 6) = Omega_(6, 7);
		Omega_(8, 0) = Omega_(0, 8); Omega_(8, 1) = Omega_(1, 8); Omega_(8, 2) = Omega_(2, 8); Omega_(8, 3) = Omega_(3, 8); Omega_(8, 4) = Omega_(4, 8); Omega_(8, 5) = Omega_(5, 8); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);

		// Q = Sum( Qi ) = Sum( [ 1, 0, -x; 0, 1, -y; -x, -y, x^2 + y^2] )
		Eigen::Matrix<double, 3, 3> Q;
		Q(0, 0) = n; Q(0, 1) = 0; Q(0, 2) = -sum_x;
		Q(1, 0) = 0; Q(1, 1) = n; Q(1, 2) = -sum_y;
		Q(2, 0) = -sum_x; Q(2, 1) = -sum_y; Q(2, 2) = sum_x2_plus_y2;

		// Qinv = inv( Q ) = inv( Sum( Qi) )
		double inv_n = 1.0 / n;
		double detQ = n * (n*sum_x2_plus_y2 - sum_y * sum_y - sum_x * sum_x);
		double point_coordinate_variance = detQ * inv_n * inv_n * inv_n;
		if (point_coordinate_variance < parameters_.point_variance_threshold)
		{
			flag_valid_ = false;
			return;
		}

		Eigen::Matrix<double, 3, 3> Qinv;
		InvertSymmetric3x3(Q, Qinv);

		// Compute P = -inv( Sum(Qi) ) * Sum( Qi*Ai ) = -Qinv * QA
		P_ = -Qinv * QA;
		// Complete Omega (i.e., Omega = Sum(A'*Qi*A') + Sum(Qi*Ai)'*P = Sum(A'*Qi*A') + Sum(Qi*Ai)'*inv(Sum(Qi))*Sum( Qi*Ai) 
		Omega_ += QA.transpose()*P_;

		// Finally, decompose Omega
		Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(Omega_, Eigen::ComputeFullU);

		U_ = svd.matrixU();
		s_ = svd.singularValues();

		// Find dimension of null space
		while (s_[7 - num_null_vectors_] < _parameters.rank_tolerance)
		{
			num_null_vectors_++;
		}

		// Dimension of null space of Omega must be <= 6
		if (++num_null_vectors_ > 6)
		{
			flag_valid_ = false;
		}

		// Point mean 
		point_mean_ << sum_X * inv_n, sum_Y*inv_n, sum_Z*inv_n;

		// Assign nearest rotation method
		if (parameters_.nearest_rotation_method == NearestRotationMethod::FOAM)
		{
			NearestRotationMatrix = NearestRotationMatrix_FOAM;
		}
		else // if ( parameters_.nearest_rotation_method == NearestRotationMethod::SVD )
		{
			NearestRotationMatrix = NearestRotationMatrix_SVD;
		}

		//以上是构造函数中的内容

		this->Solve();


		//返回旋转矩阵
		SQPSolution  sol = *this->SolutionPtr(0);

		R(0, 0) = sol.r(0);
		R(1, 0) = sol.r(1);
		R(2, 0) = sol.r(2);

		R(0, 1) = sol.r(3);
		R(1, 1) = sol.r(4);
		R(2, 1) = sol.r(5);

		R(0, 2) = sol.r(6);
		R(1, 2) = sol.r(7);
		R(2, 2) = sol.r(8);

		//平移向量
		t = sol.t;

	}


	void SqPnPSolver::ComputeRT(
		Eigen::Matrix3Xd XXw,
		Eigen::Matrix2Xd xxn,
		Eigen::Matrix3d &R,
		Eigen::Vector3d &t,
		const SolverParameters& _parameters /*= SolverParameters()*/)
	{

		std::size_t m = XXw.cols();

		std::size_t n2 = xxn.cols();


		if (m != n2)
		{
			std::cout << "error" << std::endl;
			return;
		}

		std::vector<_Point>_3dpoints;
		std::vector<_Projection> _projections;

		for (std::size_t i = 0; i < n2; i++)
		{
			Eigen::Vector3d vec3 = XXw.col(i);

			_3dpoints.push_back(_Point(vec3(0), vec3(1), vec3(2)));

			Eigen::Vector2d vec2 = xxn.col(i);
			_projections.push_back(_Projection(vec2(0), vec2(1)));

		}


		//以下是构造函数中的内容

		/*
		std::vector<_Point>_3dpoints = XXw;
		std::vector<_Projection> _projections = xxn;

		if (_3dpoints.size() != _projections.size() ||
			_3dpoints.size() < 3 || _projections.size() < 3)
		{
			flag_valid_ = false;
			return;
		}
*/

		flag_valid_ = true;
		num_null_vectors_ = -1; // set to -1 in case we never make it to the decomposition of Omega
		Omega_ = Eigen::Matrix<double, 9, 9>::Zero();
		const size_t n = _3dpoints.size();
		double sum_x = 0,
			sum_y = 0,
			sum_x2_plus_y2 = 0;

		double sum_X = 0,
			sum_Y = 0,
			sum_Z = 0;

		Eigen::Matrix<double, 3, 9> QA = Eigen::Matrix<double, 3, 9>::Zero();  // Sum( Qi*Ai )

		for (size_t i = 0; i < n; i++)
		{
			points_.push_back(_3dpoints.at(i));
			projections_.push_back(_projections.at(i));

			double x = projections_.rbegin()->vector[0],
				y = projections_.rbegin()->vector[1],
				sq_norm_m = projections_.rbegin()->vector.squaredNorm();
			sum_x += x;
			sum_y += y;
			sum_x2_plus_y2 += sq_norm_m;

			double X = points_.rbegin()->vector[0],
				Y = points_.rbegin()->vector[1],
				Z = points_.rbegin()->vector[2];
			sum_X += X;
			sum_Y += Y;
			sum_Z += Z;

			// Accumulate Omega by kronecker( Qi, Mi*Mi' ) = Q + A'*Qi*Ai. NOTE: Skipping block (3:5, 3:5) because its same as (0:2, 0:2)
			double X2 = X * X, XY = X * Y, XZ = X * Z, Y2 = Y * Y, YZ = Y * Z, Z2 = Z * Z;
			// a. Block (0:2, 0:2) populated by Mi*Mi'. NOTE: Only upper triangle
			Omega_(0, 0) += X2;
			Omega_(0, 1) += XY;
			Omega_(0, 2) += XZ;
			Omega_(1, 1) += Y2;
			Omega_(1, 2) += YZ;
			Omega_(2, 2) += Z2;

			// b. Block (0:2, 6:8) populated by -x*Mi*Mi'. NOTE: Only upper triangle
			Omega_(0, 6) += -x * X2; Omega_(0, 7) += -x * XY; Omega_(0, 8) += -x * XZ;
			Omega_(1, 7) += -x * Y2; Omega_(1, 8) += -x * YZ;
			Omega_(2, 8) += -x * Z2;
			// c. Block (3:5, 6:8) populated by -y*Mi*Mi'. NOTE: Only upper triangle
			Omega_(3, 6) += -y * X2; Omega_(3, 7) += -y * XY; Omega_(3, 8) += -y * XZ;
			Omega_(4, 7) += -y * Y2; Omega_(4, 8) += -y * YZ;
			Omega_(5, 8) += -y * Z2;
			// d. Block (6:8, 6:8) populated by (x^2+y^2)*Mi*Mi'. NOTE: Only upper triangle
			Omega_(6, 6) += sq_norm_m * X2; Omega_(6, 7) += sq_norm_m * XY; Omega_(6, 8) += sq_norm_m * XZ;
			Omega_(7, 7) += sq_norm_m * Y2; Omega_(7, 8) += sq_norm_m * YZ;
			Omega_(8, 8) += sq_norm_m * Z2;

			// Accumulating Qi*Ai in QA
			QA(0, 0) += X; QA(0, 1) += Y; QA(0, 2) += Z; 	QA(0, 6) += -x * X; QA(0, 7) += -x * Y; QA(0, 8) += -x * Z;
			QA(1, 3) += X; QA(1, 4) += Y; QA(1, 5) += Z; 	QA(1, 6) += -y * X; QA(1, 7) += -y * Y; QA(1, 8) += -y * Z;

			QA(2, 0) += -x * X; QA(2, 1) += -x * Y; QA(2, 2) += -x * Z; 	QA(2, 3) += -y * X; QA(2, 4) += -y * Y; QA(2, 5) += -y * Z;
			QA(2, 6) += sq_norm_m * X; QA(2, 7) += sq_norm_m * Y; QA(2, 8) += sq_norm_m * Z;

		}

		// Fill-in lower triangles of off-diagonal blocks (0:2, 6:8), (3:5, 6:8) and (6:8, 6:8)
		Omega_(1, 6) = Omega_(0, 7); Omega_(2, 6) = Omega_(0, 8); Omega_(2, 7) = Omega_(1, 8);
		Omega_(4, 6) = Omega_(3, 7); Omega_(5, 6) = Omega_(3, 8); Omega_(5, 7) = Omega_(4, 8);
		Omega_(7, 6) = Omega_(6, 7); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);

		// Fill-in upper triangle of block (3:5, 3:5)
		Omega_(3, 3) = Omega_(0, 0); Omega_(3, 4) = Omega_(0, 1); Omega_(3, 5) = Omega_(0, 2);
		Omega_(4, 4) = Omega_(1, 1); Omega_(4, 5) = Omega_(1, 2);
		Omega_(5, 5) = Omega_(2, 2);
		// Fill lower triangle of Omega
		Omega_(1, 0) = Omega_(0, 1);
		Omega_(2, 0) = Omega_(0, 2); Omega_(2, 1) = Omega_(1, 2);
		Omega_(3, 0) = Omega_(0, 3); Omega_(3, 1) = Omega_(1, 3); Omega_(3, 2) = Omega_(2, 3);
		Omega_(4, 0) = Omega_(0, 4); Omega_(4, 1) = Omega_(1, 4); Omega_(4, 2) = Omega_(2, 4); Omega_(4, 3) = Omega_(3, 4);
		Omega_(5, 0) = Omega_(0, 5); Omega_(5, 1) = Omega_(1, 5); Omega_(5, 2) = Omega_(2, 5); Omega_(5, 3) = Omega_(3, 5); Omega_(5, 4) = Omega_(4, 5);
		Omega_(6, 0) = Omega_(0, 6); Omega_(6, 1) = Omega_(1, 6); Omega_(6, 2) = Omega_(2, 6); Omega_(6, 3) = Omega_(3, 6); Omega_(6, 4) = Omega_(4, 6); Omega_(6, 5) = Omega_(5, 6);
		Omega_(7, 0) = Omega_(0, 7); Omega_(7, 1) = Omega_(1, 7); Omega_(7, 2) = Omega_(2, 7); Omega_(7, 3) = Omega_(3, 7); Omega_(7, 4) = Omega_(4, 7); Omega_(7, 5) = Omega_(5, 7); Omega_(7, 6) = Omega_(6, 7);
		Omega_(8, 0) = Omega_(0, 8); Omega_(8, 1) = Omega_(1, 8); Omega_(8, 2) = Omega_(2, 8); Omega_(8, 3) = Omega_(3, 8); Omega_(8, 4) = Omega_(4, 8); Omega_(8, 5) = Omega_(5, 8); Omega_(8, 6) = Omega_(6, 8); Omega_(8, 7) = Omega_(7, 8);

		// Q = Sum( Qi ) = Sum( [ 1, 0, -x; 0, 1, -y; -x, -y, x^2 + y^2] )
		Eigen::Matrix<double, 3, 3> Q;
		Q(0, 0) = n; Q(0, 1) = 0; Q(0, 2) = -sum_x;
		Q(1, 0) = 0; Q(1, 1) = n; Q(1, 2) = -sum_y;
		Q(2, 0) = -sum_x; Q(2, 1) = -sum_y; Q(2, 2) = sum_x2_plus_y2;

		// Qinv = inv( Q ) = inv( Sum( Qi) )
		double inv_n = 1.0 / n;
		double detQ = n * (n*sum_x2_plus_y2 - sum_y * sum_y - sum_x * sum_x);
		double point_coordinate_variance = detQ * inv_n * inv_n * inv_n;
		if (point_coordinate_variance < parameters_.point_variance_threshold)
		{
			flag_valid_ = false;
			return;
		}

		Eigen::Matrix<double, 3, 3> Qinv;
		InvertSymmetric3x3(Q, Qinv);

		// Compute P = -inv( Sum(Qi) ) * Sum( Qi*Ai ) = -Qinv * QA
		P_ = -Qinv * QA;
		// Complete Omega (i.e., Omega = Sum(A'*Qi*A') + Sum(Qi*Ai)'*P = Sum(A'*Qi*A') + Sum(Qi*Ai)'*inv(Sum(Qi))*Sum( Qi*Ai) 
		Omega_ += QA.transpose()*P_;

		// Finally, decompose Omega
		Eigen::JacobiSVD<Eigen::Matrix<double, 9, 9>> svd(Omega_, Eigen::ComputeFullU);

		U_ = svd.matrixU();
		s_ = svd.singularValues();

		// Find dimension of null space
		while (s_[7 - num_null_vectors_] < _parameters.rank_tolerance)
		{
			num_null_vectors_++;
		}

		// Dimension of null space of Omega must be <= 6
		if (++num_null_vectors_ > 6)
		{
			flag_valid_ = false;
		}

		// Point mean 
		point_mean_ << sum_X * inv_n, sum_Y*inv_n, sum_Z*inv_n;

		// Assign nearest rotation method
		if (parameters_.nearest_rotation_method == NearestRotationMethod::FOAM)
		{
			NearestRotationMatrix = NearestRotationMatrix_FOAM;
		}
		else // if ( parameters_.nearest_rotation_method == NearestRotationMethod::SVD )
		{
			NearestRotationMatrix = NearestRotationMatrix_SVD;
		}

		//以上是构造函数中的内容

		this->Solve();


		//返回旋转矩阵
		SQPSolution  sol = *this->SolutionPtr(0);

		R(0, 0) = sol.r(0);
		R(1, 0) = sol.r(1);
		R(2, 0) = sol.r(2);

		R(0, 1) = sol.r(3);
		R(1, 1) = sol.r(4);
		R(2, 1) = sol.r(5);

		R(0, 2) = sol.r(6);
		R(1, 2) = sol.r(7);
		R(2, 2) = sol.r(8);

		//平移向量
		t = sol.t;






	}




	//
	// Solve the PnP 
	bool SqPnPSolver::Solve()
	{
		if (!flag_valid_)
		{
			return false;
		}

		double min_sq_error = std::numeric_limits<double>::max();
		int num_eigen_points = num_null_vectors_ > 0 ? num_null_vectors_ : 1;

		// clear solutions
		num_solutions_ = 0;

		for (int i = 9 - num_eigen_points; i < 9; i++)
		{
			// NOTE: No need to scale by sqrt(3) here, but better be there for other computations (i.e., orthogonality test)
			const Eigen::Matrix<double, 9, 1> e = SQRT3 * Eigen::Map<Eigen::Matrix<double, 9, 1>>(U_.block<9, 1>(0, i).data());
			double orthogonality_sq_error = OrthogonalityError(e);
			// Find nearest rotation vector
			SQPSolution solution[2];

			// Avoid SQP if e is orthogonal
			if (orthogonality_sq_error < parameters_.orthogonality_squared_error_threshold)
			{
				solution[0].r_hat = Determinant9x1(e) * e;
				solution[0].t = P_ * solution[0].r_hat;
				solution[0].num_iterations = 0;

				HandleSolution(solution[0], min_sq_error);
			}
			else
			{
				NearestRotationMatrix(e, solution[0].r);
				solution[0] = RunSQP(solution[0].r);
				solution[0].t = P_ * solution[0].r_hat;
				HandleSolution(solution[0], min_sq_error);

				NearestRotationMatrix(-e, solution[1].r);
				solution[1] = RunSQP(solution[1].r);
				solution[1].t = P_ * solution[1].r_hat;
				HandleSolution(solution[1], min_sq_error);
			}
		}

		int c = 1;
		while (min_sq_error > 3 * s_[9 - num_eigen_points - c] && 9 - num_eigen_points - c > 0)
		{
			int index = 9 - num_eigen_points - c;

			const Eigen::Matrix<double, 9, 1> e = Eigen::Map<Eigen::Matrix<double, 9, 1>>(U_.block<9, 1>(0, index).data());
			SQPSolution solution[2];

			NearestRotationMatrix(e, solution[0].r);
			solution[0] = RunSQP(solution[0].r);
			solution[0].t = P_ * solution[0].r_hat;
			HandleSolution(solution[0], min_sq_error);

			NearestRotationMatrix(-e, solution[1].r);
			solution[1] = RunSQP(solution[1].r);
			solution[1].t = P_ * solution[1].r_hat;
			HandleSolution(solution[1], min_sq_error);

			c++;
		}

		return true;
	}


	//
	// Run sequential quadratic programming on orthogonal matrices
	SQPSolution SqPnPSolver::RunSQP(const Eigen::Matrix<double, 9, 1>& r0)
	{
		Eigen::Matrix<double, 9, 1> r = r0;

		double delta_squared_norm = std::numeric_limits<double>::max();
		Eigen::Matrix<double, 9, 1> delta;
		int step = 0;

		while (delta_squared_norm > parameters_.sqp_squared_tolerance && step++ < parameters_.sqp_max_iteration)
		{
			SolveSQPSystem(r, delta);
			r += delta;
			delta_squared_norm = delta.squaredNorm();
		}

		SQPSolution solution;
		solution.num_iterations = step;
		solution.r = r;
		// clear the estimate and/or flip the matrix sign if necessary
		double det_r = Determinant9x1(solution.r);
		if (det_r < 0)
		{
			solution.r = -r;
			det_r = -det_r;
		}

		if (det_r > parameters_.sqp_det_threshold)
		{
			NearestRotationMatrix(solution.r, solution.r_hat);
		}
		else
		{
			solution.r_hat = solution.r;
		}


		return solution;
	}

	//
	// Solve the SQP system efficiently
	void SqPnPSolver::SolveSQPSystem(const Eigen::Matrix<double, 9, 1>& r, Eigen::Matrix<double, 9, 1>& delta)
	{
		double sqnorm_r1 = r[0] * r[0] + r[1] * r[1] + r[2] * r[2],
			sqnorm_r2 = r[3] * r[3] + r[4] * r[4] + r[5] * r[5],
			sqnorm_r3 = r[6] * r[6] + r[7] * r[7] + r[8] * r[8];

		double dot_r1r2 = r[0] * r[3] + r[1] * r[4] + r[2] * r[5],
			dot_r1r3 = r[0] * r[6] + r[1] * r[7] + r[2] * r[8],
			dot_r2r3 = r[3] * r[6] + r[4] * r[7] + r[5] * r[8];

		// Obtain 6D normal (H) and 3D null space of the constraint Jacobian-J at the estimate (r)
		// NOTE: Thsi is done via Gram-Schmidt orthogoalization
		Eigen::Matrix<double, 9, 3> N;  // Null space of J
		Eigen::Matrix<double, 9, 6> H;  // Row space of J
		Eigen::Matrix<double, 6, 6> JH; // The lower triangular matrix J*Q

		RowAndNullSpace(r, H, N, JH);

		// Great, now if delta = H*x + N*y, we first compute x by solving:
		// 
		//              (J*H)*x = g
		//
		// where g is the constraint vector g = [   1 - norm(r1)^2;
		// 					     	   1 - norm(r2)^2;
		//					     	   1 - norm(r3)^2;
		//					           -r1'*r2; 
		//						   -r2'*r3; 
		//						   -r1'*r3 ];
		Eigen::Matrix<double, 6, 1> g;
		g[0] = 1 - sqnorm_r1; g[1] = 1 - sqnorm_r2; g[2] = 1 - sqnorm_r3; g[3] = -dot_r1r2; g[4] = -dot_r2r3; g[5] = -dot_r1r3;

		Eigen::Matrix<double, 6, 1> x;
		x[0] = g[0] / JH(0, 0);
		x[1] = g[1] / JH(1, 1);
		x[2] = g[2] / JH(2, 2);
		x[3] = (g[3] - JH(3, 0)*x[0] - JH(3, 1)*x[1]) / JH(3, 3);
		x[4] = (g[4] - JH(4, 1)*x[1] - JH(4, 2)*x[2] - JH(4, 3)*x[3]) / JH(4, 4);
		x[5] = (g[5] - JH(5, 0)*x[0] - JH(5, 2)*x[2] - JH(5, 3)*x[3] - JH(5, 4)*x[4]) / JH(5, 5);

		// Now obtain the component of delta in the row space of E as delta_h = Q'*x and assign straint into delta
		delta = H * x;

		// Finally, solve for y from W*y = ksi , where matrix W and vector ksi are :
		//
		// W = N'*Omega*N and ksi = -N'*Omega*( r + delta_h );
		Eigen::Matrix<double, 3, 9> NtOmega = N.transpose() * Omega_;
		Eigen::Matrix<double, 3, 3> W = NtOmega * N, Winv;
		InvertSymmetric3x3(W, Winv); // NOTE: This maybe also analytical with Eigen, but hey...

		Eigen::Matrix<double, 3, 1> y = -Winv * NtOmega * (delta + r);

		// FINALLY, accumulate delta with component in tangent space (delta_n)
		delta += N * y;
	}


	//
	// Compute the 3D null space (N) and 6D normal space (H) of the constraint Jacobian at a 9D vector r 
	// (r is not necessarilly a rotation but it must represent an rank-3 matrix )
	// NOTE: K is lower-triangular, so upper triangle may contain trash (is not filled by the function)...
	void SqPnPSolver::RowAndNullSpace(
		const Eigen::Matrix<double, 9, 1>& r,
		Eigen::Matrix<double, 9, 6>& H, // Row space 
		Eigen::Matrix<double, 9, 3>& N, // Null space
		Eigen::Matrix<double, 6, 6>& K,  // J*Q (J - Jacobian of constraints)
		const double& norm_threshold // Used to discard columns of Pn when finding null space
	) // threshold for column vector norm (of Pn)
	{
		// Applying Gram-Schmidt orthogonalization on the Jacobian. 
		// The steps are fixed here to take advantage of the sparse form of the matrix
		//
		H = Eigen::Matrix<double, 9, 6>::Zero();

		// 1. q1
		double norm_r1 = sqrt(r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
		double inv_norm_r1 = norm_r1 > 1e-5 ? 1.0 / norm_r1 : 0.0;
		H(0, 0) = r[0] * inv_norm_r1; H(1, 0) = r[1] * inv_norm_r1; H(2, 0) = r[2] * inv_norm_r1;
		K(0, 0) = 2 * norm_r1;

		// 2. q2 
		double norm_r2 = sqrt(r[3] * r[3] + r[4] * r[4] + r[5] * r[5]);
		double inv_norm_r2 = 1.0 / norm_r2;
		H(3, 1) = r[3] * inv_norm_r2; H(4, 1) = r[4] * inv_norm_r2; H(5, 1) = r[5] * inv_norm_r2;
		K(1, 0) = 0; K(1, 1) = 2 * norm_r2;

		// 3. q3 = (r3'*q2)*q2 - (r3'*q1)*q1 ; q3 = q3/norm(q3)
		double norm_r3 = sqrt(r[6] * r[6] + r[7] * r[7] + r[8] * r[8]);
		double inv_norm_r3 = 1.0 / norm_r3;
		H(6, 2) = r[6] * inv_norm_r3; H(7, 2) = r[7] * inv_norm_r3; H(8, 2) = r[8] * inv_norm_r3;
		K(2, 0) = K(2, 1) = 0; K(2, 2) = 2 * norm_r3;

		// 4. q4
		double dot_j4q1 = r[3] * H(0, 0) + r[4] * H(1, 0) + r[5] * H(2, 0),
			dot_j4q2 = r[0] * H(3, 1) + r[1] * H(4, 1) + r[2] * H(5, 1);

		H(0, 3) = r[3] - dot_j4q1 * H(0, 0); H(1, 3) = r[4] - dot_j4q1 * H(1, 0); H(2, 3) = r[5] - dot_j4q1 * H(2, 0);
		H(3, 3) = r[0] - dot_j4q2 * H(3, 1); H(4, 3) = r[1] - dot_j4q2 * H(4, 1); H(5, 3) = r[2] - dot_j4q2 * H(5, 1);
		double inv_norm_j4 = 1.0 / sqrt(H(0, 3)*H(0, 3) + H(1, 3)*H(1, 3) + H(2, 3)*H(2, 3) +
			H(3, 3)*H(3, 3) + H(4, 3)*H(4, 3) + H(5, 3)*H(5, 3));

		H(0, 3) *= inv_norm_j4; H(1, 3) *= inv_norm_j4; H(2, 3) *= inv_norm_j4;
		H(3, 3) *= inv_norm_j4; H(4, 3) *= inv_norm_j4; H(5, 3) *= inv_norm_j4;

		K(3, 0) = r[3] * H(0, 0) + r[4] * H(1, 0) + r[5] * H(2, 0); K(3, 1) = r[0] * H(3, 1) + r[1] * H(4, 1) + r[2] * H(5, 1);
		K(3, 2) = 0; K(3, 3) = r[3] * H(0, 3) + r[4] * H(1, 3) + r[5] * H(2, 3) + r[0] * H(3, 3) + r[1] * H(4, 3) + r[2] * H(5, 3);

		// 5. q5
		double dot_j5q2 = r[6] * H(3, 1) + r[7] * H(4, 1) + r[8] * H(5, 1),
			dot_j5q3 = r[3] * H(6, 2) + r[4] * H(7, 2) + r[5] * H(8, 2),
			dot_j5q4 = r[6] * H(3, 3) + r[7] * H(4, 3) + r[8] * H(5, 3);

		H(0, 4) = -dot_j5q4 * H(0, 3); 			  H(1, 4) = -dot_j5q4 * H(1, 3); 			H(2, 4) = -dot_j5q4 * H(2, 3);
		H(3, 4) = r[6] - dot_j5q2 * H(3, 1) - dot_j5q4 * H(3, 3); H(4, 4) = r[7] - dot_j5q2 * H(4, 1) - dot_j5q4 * H(4, 3); H(5, 4) = r[8] - dot_j5q2 * H(5, 1) - dot_j5q4 * H(5, 3);
		H(6, 4) = r[3] - dot_j5q3 * H(6, 2); H(7, 4) = r[4] - dot_j5q3 * H(7, 2); H(8, 4) = r[5] - dot_j5q3 * H(8, 2);

		H.block<9, 1>(0, 4) /= H.col(4).norm();

		K(4, 0) = 0; K(4, 1) = r[6] * H(3, 1) + r[7] * H(4, 1) + r[8] * H(5, 1); K(4, 2) = r[3] * H(6, 2) + r[4] * H(7, 2) + r[5] * H(8, 2);
		K(4, 3) = r[6] * H(3, 3) + r[7] * H(4, 3) + r[8] * H(5, 3);
		K(4, 4) = r[6] * H(3, 4) + r[7] * H(4, 4) + r[8] * H(5, 4) + r[3] * H(6, 4) + r[4] * H(7, 4) + r[5] * H(8, 4);


		// 4. q6
		double dot_j6q1 = r[6] * H(0, 0) + r[7] * H(1, 0) + r[8] * H(2, 0),
			dot_j6q3 = r[0] * H(6, 2) + r[1] * H(7, 2) + r[2] * H(8, 2),
			dot_j6q4 = r[6] * H(0, 3) + r[7] * H(1, 3) + r[8] * H(2, 3),
			dot_j6q5 = r[0] * H(6, 4) + r[1] * H(7, 4) + r[2] * H(8, 4) + r[6] * H(0, 4) + r[7] * H(1, 4) + r[8] * H(2, 4);

		H(0, 5) = r[6] - dot_j6q1 * H(0, 0) - dot_j6q4 * H(0, 3) - dot_j6q5 * H(0, 4);
		H(1, 5) = r[7] - dot_j6q1 * H(1, 0) - dot_j6q4 * H(1, 3) - dot_j6q5 * H(1, 4);
		H(2, 5) = r[8] - dot_j6q1 * H(2, 0) - dot_j6q4 * H(2, 3) - dot_j6q5 * H(2, 4);

		H(3, 5) = -dot_j6q5 * H(3, 4) - dot_j6q4 * H(3, 3);
		H(4, 5) = -dot_j6q5 * H(4, 4) - dot_j6q4 * H(4, 3);
		H(5, 5) = -dot_j6q5 * H(5, 4) - dot_j6q4 * H(5, 3);

		H(6, 5) = r[0] - dot_j6q3 * H(6, 2) - dot_j6q5 * H(6, 4);
		H(7, 5) = r[1] - dot_j6q3 * H(7, 2) - dot_j6q5 * H(7, 4);
		H(8, 5) = r[2] - dot_j6q3 * H(8, 2) - dot_j6q5 * H(8, 4);

		H.block<9, 1>(0, 5) /= H.col(5).norm();

		K(5, 0) = r[6] * H(0, 0) + r[7] * H(1, 0) + r[8] * H(2, 0); K(5, 1) = 0; K(5, 2) = r[0] * H(6, 2) + r[1] * H(7, 2) + r[2] * H(8, 2);
		K(5, 3) = r[6] * H(0, 3) + r[7] * H(1, 3) + r[8] * H(2, 3); K(5, 4) = r[6] * H(0, 4) + r[7] * H(1, 4) + r[8] * H(2, 4) + r[0] * H(6, 4) + r[1] * H(7, 4) + r[2] * H(8, 4);
		K(5, 5) = r[6] * H(0, 5) + r[7] * H(1, 5) + r[8] * H(2, 5) + r[0] * H(6, 5) + r[1] * H(7, 5) + r[2] * H(8, 5);

		// Great! Now H is an orthogonalized, sparse basis of the Jacobian row space and K is filled.
		//
		// Now get a projector onto the null space of H:
		const Eigen::Matrix<double, 9, 9> Pn = Eigen::Matrix<double, 9, 9>::Identity() - (H*H.transpose());

		// Now we need to pick 3 columns of P with non-zero norm (> 0.3) and some angle between them (> 0.3).
		//
		// Find the 3 columns of Pn with largest norms
		int index1 = -1,
			index2 = -1,
			index3 = -1;
		double  max_norm1 = std::numeric_limits<double>::min(),
			min_dot12 = std::numeric_limits<double>::max(),
			min_dot1323 = std::numeric_limits<double>::max();


		double col_norms[9];
		for (int i = 0; i < 9; i++)
		{
			col_norms[i] = Pn.col(i).norm();
			if (col_norms[i] >= norm_threshold)
			{
				if (max_norm1 < col_norms[i])
				{
					max_norm1 = col_norms[i];
					index1 = i;
				}
			}
		}

		const auto& v1 = Pn.block<9, 1>(0, index1);
		N.block<9, 1>(0, 0) = v1 * (1.0 / max_norm1);

		for (int i = 0; i < 9; i++)
		{
			if (i == index1)
			{
				continue;
			}

			if (col_norms[i] >= norm_threshold)
			{
				double cos_v1_x_col = std::fabs(Pn.col(i).dot(v1) / col_norms[i]);

				if (cos_v1_x_col <= min_dot12)
				{
					index2 = i;
					min_dot12 = cos_v1_x_col;
				}
			}
		}

		const auto& v2 = Pn.block<9, 1>(0, index2);
		N.block<9, 1>(0, 1) = v2 - v2.dot(N.col(0)) * N.col(0);
		N.block<9, 1>(0, 1) /= N.col(1).norm();

		for (int i = 0; i < 9; i++)
		{
			if (i == index2 || i == index1) continue;
			if (col_norms[i] >= norm_threshold)
			{
				double cos_v1_x_col = std::fabs(Pn.col(i).dot(v1) / col_norms[i]);
				double cos_v2_x_col = std::fabs(Pn.col(i).dot(v2) / col_norms[i]);

				if (cos_v1_x_col + cos_v2_x_col <= min_dot1323)
				{
					index3 = i;
					min_dot1323 = cos_v2_x_col + cos_v2_x_col;
				}
			}
		}

		// Now orthogonalize the remaining 2 vectors v2, v3 into N
		const auto& v3 = Pn.block<9, 1>(0, index3);

		N.block<9, 1>(0, 2) = v3 - (v3.dot(N.col(1)) * N.col(1)) - (v3.dot(N.col(0)) * N.col(0));
		N.block<9, 1>(0, 2) /= N.col(2).norm();

	}


	std::size_t str2Int(std::string str)
	{
		if (str.empty())
		{
			std::cout << "str is empty !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;
		}

		std::istringstream iss(str);

		std::size_t cnt;

		iss >> cnt;

		return cnt;

	}



	double str2Double(std::string str)
	{
		if (str.empty())
		{
			std::cout << "str is empty !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;
			return -1;
		}

		std::istringstream iss(str);

		double cnt;

		iss >> cnt;

		return cnt;

	}


	void split(std::string lineStr, std::vector<std::string> &vecStr)
	{

		if (lineStr.empty())
		{
			std::cout << "lineStr is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;
			return;

		}

		lineStr = boost::trim_copy(lineStr);

		boost::split(vecStr, lineStr, boost::is_any_of("\t"), boost::token_compress_off);

	}



	void read3D(std::string fileName, std::vector<_Point> &XXw)
	{
		if (fileName.empty() == true)
		{
			std::cout << "FileName is empty !, in line: "
				<< __LINE__ << ", at file: " << __FILE__ << std::endl;

			return;
		}

		std::ifstream file;
		file.open(fileName, std::ifstream::in);

		if (!file.is_open())
		{
			std::cout << "Failed to open file !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return;

		}

		std::string numStr;
		std::getline(file, numStr);

		std::size_t cnt = str2Int(numStr);


		std::string lineStr;

		std::vector<std::string >vecStr;

		std::size_t counter = 0;

		while (std::getline(file, lineStr) && counter < cnt)
		{
			//分割字符串
			split(lineStr, vecStr);

			double e1 = str2Double(vecStr[0]);
			double e2 = str2Double(vecStr[1]);
			double e3 = str2Double(vecStr[2]);

			XXw.push_back(cvg::_Point(e1, e2, e3));

			//
			counter++;

		}

		//关闭文件
		file.close();


	}//end of read2D;


	void read2D(std::string fileName, std::vector<_Projection> &xxn)
	{

		if (fileName.empty())
		{

			std::cout << "fileName is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return;
		}


		//打开文件
		std::ifstream file;
		file.open(fileName, std::ifstream::in);

		if (!file.is_open())
		{
			std::cout << "Failed to open file !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return;
		}

		std::string numStr;

		std::getline(file, numStr);

		std::size_t cnt = str2Int(numStr);


		std::string lineStr;
		std::vector<std::string >vecStr;
		std::size_t counter = 0;

		while (std::getline(file, lineStr) && counter < cnt)
		{
			//分割字符串
			split(lineStr, vecStr);

			double e1 = str2Double(vecStr[0]);
			double e2 = str2Double(vecStr[1]);

			xxn.push_back(cvg::_Projection(e1, e2));

			counter++;
		}

		//关闭文件
		file.close();

	}//end of read2D;


	void readData(
		std::string file3D,
		std::string file2D,
		std::vector<_Point> &XXw,
		std::vector<_Projection> &xxn)
	{

		//清空缓冲区
		XXw.clear();
		xxn.clear();

		//读取3D,2D点信息
		read3D(file3D, XXw);

		read2D(file2D, xxn);

	}

}
