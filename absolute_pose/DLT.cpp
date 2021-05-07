

#pragma warning(disable:4996)


#include "DLT.hpp"

#include <iostream>


namespace cvg
{


	CDLT::CDLT()
	{

	}


	CDLT::~CDLT()
	{

	}

	int CDLT::ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
		Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, double tol /*= 1e-8*/)
	{
		if (MatXYZ.cols() == 0 || MatXYZ.rows() == 0)
		{
			std::cout << "MatXYZ is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		if (Matxy.cols() == 0 || Matxy.rows() == 0)
		{
			std::cout << "Matxy is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}

		//获取点数.
		std::size_t num = Matxy.cols();

		Eigen::MatrixXd matD = Eigen::MatrixXd::Zero(num * 2, 12);


		double xi = 0.0, yi = 0.0, zi = 0.0;
		double ui = 0.0, vi = 0.0;


		for (std::size_t i = 0; i < num; i++)
		{
			//获取三维点坐标.
			xi = MatXYZ.col(i)(0);
			yi = MatXYZ.col(i)(1);
			zi = MatXYZ.col(i)(2);

			//获取二维点坐标.
			ui = Matxy.col(i)(0);
			vi = Matxy.col(i)(1);

			//矩阵第2*i行.
			matD.row(2 * i)(0) = xi;
			matD.row(2 * i)(1) = yi;
			matD.row(2 * i)(2) = zi;

			matD.row(2 * i)(6) = -ui * xi;
			matD.row(2 * i)(7) = -ui * yi;
			matD.row(2 * i)(8) = -ui * zi;

			matD.row(2 * i)(9) = 1;
			matD.row(2 * i)(11) = -ui;

			//矩阵第2*i+1行.
			matD.row(2 * i + 1)(3) = xi;
			matD.row(2 * i + 1)(4) = yi;
			matD.row(2 * i + 1)(5) = zi;

			matD.row(2 * i + 1)(6) = -vi * xi;
			matD.row(2 * i + 1)(7) = -vi * yi;
			matD.row(2 * i + 1)(8) = -vi * zi;

			matD.row(2 * i + 1)(10) = 1;
			matD.row(2 * i + 1)(11) = -vi;

		}

		//std::cout << "xxxx" << std::endl;

		//
		Eigen::MatrixXd matDD = matD.transpose() * matD;

		Eigen::MatrixXd matV;
		EigenSolver(matDD, matV);



		Eigen::VectorXd vecV = matV.col(0);

		vecV = vecV / vecV.segment<3>(6).norm();

		vecV = vecV * sign(vecV(11));

		mR(0, 0) = vecV(0);
		mR(0, 1) = vecV(1);
		mR(0, 2) = vecV(2);

		mR(1, 0) = vecV(3);
		mR(1, 1) = vecV(4);
		mR(1, 2) = vecV(5);

		mR(2, 0) = vecV(6);
		mR(2, 1) = vecV(7);
		mR(2, 2) = vecV(8);

		//mT = vecV.segment<9>(3);
		mT = vecV.tail<3>();


		Eigen::MatrixXd matTemp = Eigen::MatrixXd::Zero(3, num);

		for (std::size_t i = 0; i < num; i++)
		{
			matTemp.col(i) = mT;
		}


		Eigen::Matrix3Xd matXXc = mR * MatXYZ + matTemp;


		CalcPose(mR, mT, MatXYZ, matXXc);




		return 1;

	}

	int CDLT::CalcPose(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
		Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix3Xd &Matxy)
	{

		if (MatXYZ.cols() == 0 || MatXYZ.rows() == 0)
		{
			std::cout << "MatXYZ is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}

		if (Matxy.cols() == 0 || Matxy.rows() == 0)
		{
			std::cout << "Matxy is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;
			return -1;
		}


		std::size_t n = Matxy.cols();

		Eigen::Matrix3Xd matX = MatXYZ;
		Eigen::Matrix3Xd matY = Matxy;

		Eigen::MatrixXd matK = Eigen::MatrixXd::Identity(n, n) - (Eigen::MatrixXd::Ones(n, n) / n);

		Eigen::Vector3d ux = matX.rowwise().mean();
		Eigen::Vector3d uy = matY.rowwise().mean();


		double sigmx2 = ((matX *matK).array().square()).rowwise().sum().mean();



		Eigen::Matrix3d SXY = matY * matK * matX.transpose() / n;


		Eigen::JacobiSVD<Eigen::Matrix3d> svd(SXY, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::Matrix3d matU = svd.matrixU();
		Eigen::Matrix3d matV = svd.matrixV();

		Eigen::Matrix3d matD = svd.singularValues().asDiagonal();

		Eigen::Matrix3d matS = Eigen::Matrix3d::Identity();

		if (matS.determinant() < 0)
		{
			matS(2, 2) = -1;
		}

		mR = matU * matS * matV.transpose();

		double c2 = (matD* matS).trace() / sigmx2;

		mT = uy - c2 * mR *ux;

		Eigen::Vector3d vecX = mR.col(0);
		Eigen::Vector3d vecY = mR.col(1);
		Eigen::Vector3d vecZ = mR.col(2);

		if ((vecX.cross(vecY) - vecZ).norm() > 2e-2)
		{
			mR.col(2) = -vecZ;
		}


		return 1;
	}

	int CDLT::sign(double x)
	{
		int flag = 0;

		if (x > 0)
		{
			flag = 1;
		}
		else if (x == 0)
		{
			flag = 0;
		}
		else
		{
			flag = -1;
		}


		return flag;

	}

	int CDLT::EigenSolver(Eigen::MatrixXd &mat, Eigen::MatrixXd &matV)
	{
		if (mat.cols() == 0 || mat.rows() == 0)
		{
			std::cout << "mat is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		//求解特征向量.
		Eigen::EigenSolver<Eigen::MatrixXd> eig(mat);

		matV = eig.pseudoEigenvectors();

		return 1;

	}

}//namespace cvg;
