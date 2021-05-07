
#pragma warning(disable:4996)

#include <iostream>

#include <Eigen/Eigen>


#include "HOMO.hpp"



namespace cvg
{



	CHOMO::CHOMO()
	{

	}


	CHOMO::~CHOMO()
	{

	}

	int CHOMO::ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
		Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, double tol /*= 1e-8*/)
	{


		//计算点数.
		std::size_t num = MatXYZ.cols();

		Eigen::MatrixXd matA = Eigen::MatrixXd::Zero(num * 2, 9);

		//取三维点的前两行.
		Eigen::MatrixXd matXYI = MatXYZ.topRows<2>();

		Eigen::MatrixXd matUVI = Matxy;

		double xi = 0.0, yi = 0.0;
		double ui = 0.0, vi = 0.0;




		for (std::size_t i = 0; i < num; i++)
		{

			xi = matXYI.col(i)(0);
			yi = matXYI.col(i)(1);

			ui = matUVI.col(i)(0);
			vi = matUVI.col(i)(1);


			//A矩阵的行(1).

			matA.row(2 * i)(0) = xi;
			matA.row(2 * i)(1) = yi;
			matA.row(2 * i)(2) = 1;

			matA.row(2 * i)(6) = -ui * xi;
			matA.row(2 * i)(7) = -ui * yi;
			matA.row(2 * i)(8) = -ui;

			//A矩阵的行(2).

			matA.row(2 * i + 1)(3) = xi;
			matA.row(2 * i + 1)(4) = yi;
			matA.row(2 * i + 1)(5) = 1;

			matA.row(2 * i + 1)(6) = -vi * xi;
			matA.row(2 * i + 1)(7) = -vi * yi;
			matA.row(2 * i + 1)(8) = -vi;



		}

		Eigen::MatrixXd matAtA = matA.transpose() * matA;

		//求取矩阵matV的特征向量.
		Eigen::MatrixXd matV;
		EigenSolver(matAtA, matV);

		Eigen::VectorXd vecV = matV.col(0);

		std::size_t n = vecV.size();

		vecV = vecV * sign(vecV(n - 1));



		Eigen::Matrix3d matAA;
		matAA(0, 0) = vecV(0);
		matAA(0, 1) = vecV(1);
		matAA(0, 2) = vecV(2);
		matAA(1, 0) = vecV(3);
		matAA(1, 1) = vecV(4);
		matAA(1, 2) = vecV(5);
		matAA(2, 0) = vecV(6);
		matAA(2, 1) = vecV(7);
		matAA(2, 2) = vecV(8);

		matAA = matAA / matAA.col(0).norm();

		mR = matAA;
		mT = matAA.col(2);
		mR.col(2) = matAA.col(0).cross(matAA.col(1));

		Eigen::JacobiSVD<Eigen::Matrix3d> svd(mR, Eigen::ComputeFullU | Eigen::ComputeFullV);

		mR = svd.matrixU() * svd.matrixV().transpose();



		return 1;

	}

	int CHOMO::EigenSolver(Eigen::MatrixXd &mat, Eigen::MatrixXd &matV)
	{

		if (mat.cols() == 0 || mat.rows() == 0)
		{
			return -1;
		}


		Eigen::EigenSolver<Eigen::MatrixXd> eig(mat);

		matV = eig.pseudoEigenvectors();


		return 1;


	}

	int CHOMO::sign(double x)
	{
		if (x > 0)
		{
			return 1;

		}
		else if (x == 0)
		{
			return 0;
		}
		else
		{
			return -1;
		}

	}


}//namespace cvg;
