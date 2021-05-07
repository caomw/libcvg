
#pragma warning(disable:4996)

#include <iostream>

#include <Eigen/Eigen>

#include "UP2P.h"
#include "univariate.h"

namespace cvg
{
	CUP2P::CUP2P()
	{

	}

	CUP2P::~CUP2P()
	{

	}


	void CUP2P::ComputeRT(
		Eigen::Matrix3Xd matXYZ,
		Eigen::Matrix2Xd matXY,
		Eigen::Matrix3d &R,
		Eigen::Vector3d &t)
	{
		//计算三维点的数量
		std::size_t m1 = matXYZ.rows();
		std::size_t n1 = matXYZ.cols();


		//计算二维点的数量
		std::size_t m2 = matXY.rows();
		std::size_t n2 = matXY.cols();

		if (n1 != n2)
		{
			std::cout << "input data is incorrect. in line: "
				<< __LINE__ << ", at file: " << __FILE__ << std::endl;
			return;
		}

		Eigen::Matrix3Xd matXY2 = Eigen::Matrix3Xd::Ones(3, n2);
		matXY2.topRows(2) = matXY;


		//转换三维点;
		std::vector<Eigen::Vector3d>X;

		for (std::size_t i = 0; i < m1; i++)
		{
			X.push_back(matXYZ.col(i));
		}


		//转换二维点
		std::vector<Eigen::Vector3d> x;
		for (std::size_t i = 0; i < n2; i++)
		{
			x.push_back(matXY2.col(i));
		}


		//
		CameraPoseVector pose;
		int  n = this->up2p(X, x, &pose);

		R = pose[0].R;
		t = pose[0].t;
	}


	int CUP2P::up2p(
		const std::vector<Eigen::Vector3d> &x,
		const std::vector<Eigen::Vector3d> &X,
		CameraPoseVector *output)
	{
		Eigen::Matrix<double, 4, 4> A;
		Eigen::Matrix<double, 4, 2> b;

		A << -x[0](2), 0, x[0](0), X[0](0) * x[0](2) + X[0](2) * x[0](0), 0, -x[0](2), x[0](1), X[0](1) * x[0](2) + X[0](2) * x[0](1), -x[1](2), 0, x[1](0), X[1](0) * x[1](2) + X[1](2) * x[1](0), 0, -x[1](2), x[1](1), X[1](1) * x[1](2) + X[1](2) * x[1](1);
		b << -2.0 * X[0](1) * x[0](2), X[0](0) * x[0](2) - X[0](2) * x[0](0), 2.0 * X[0](0) * x[0](2), X[0](1) * x[0](2) - X[0](2) * x[0](1), -2.0 * X[1](1) * x[1](2), X[1](0) * x[1](2) - X[1](2) * x[1](0), 2.0 * X[1](0) * x[1](2), X[1](1) * x[1](2) - X[1](2) * x[1](1);

		//b = A.partialPivLu().solve(b);
		b = A.inverse() * b;

		double c2 = -b(3, 0);
		double c3 = -b(3, 1);

		double qq[2];
		int sols = univariate::solve_quadratic_real(1.0, c2, c3, qq);

		for (int i = 0; i < sols; ++i)
		{
			CameraPose pose;

			double q = qq[i];
			double q2 = q * q;
			double cq = (1 - q2) / (1 + q2);
			double sq = 2 * q / (1 + q2);

			pose.R.setIdentity();
			pose.R(0, 0) = cq;
			pose.R(0, 1) = -sq;
			pose.R(1, 0) = sq;
			pose.R(1, 1) = cq;

			pose.t = b.block<3, 1>(0, 0) * q + b.block<3, 1>(0, 1);
			pose.t /= (1 + q2);

			output->push_back(pose);
		}

		return sols;
	}

}
