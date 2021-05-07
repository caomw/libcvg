
#pragma warning(disable:4996)

#include <iostream>

#include <Eigen/Eigen>

#include "Umeyama.h"


namespace cvg
{

	CUmeyama::CUmeyama()
	{

	}

	CUmeyama::~CUmeyama()
	{

	}


	void CUmeyama::ComputeRT(
		Eigen::Matrix3Xd matU, Eigen::Matrix2Xd matv,
		Eigen::Matrix3d &poseR, Eigen::Vector3d &poseT)
	{


		if (matU.cols() != matv.cols())
		{
			std::cout << "input data is not consistency. in line: "
				<< __LINE__ << ", at file: " << __FILE__ << std::endl;

			return;
		}

		if (matU.rows() != 3 || matv.rows() != 2)
		{
			std::cout << "input data is incorrect, in line: "
				<< __LINE__ << ", at file: " << __FILE__ << std::endl;

			return;
		}


		Eigen::Matrix3Xd Y = matU;

		std::size_t m = matU.rows();
		std::size_t n = matU.cols();

		Eigen::Matrix3Xd X = Eigen::Matrix3Xd::Ones(m, n);
		X.topRows(2) = matv;



		//////

		Eigen::Vector3d mean_X = X.rowwise().mean();
		Eigen::Vector3d mean_Y = Y.rowwise().mean();

		Eigen::Matrix3Xd X_demean = X.colwise() - mean_X;
		Eigen::Matrix3Xd Y_demean = Y.colwise() - mean_Y;

		Eigen::MatrixXd sigma = 1.0 / n * Y_demean*X_demean.transpose();

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(sigma, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::MatrixXd U = svd.matrixU();
		Eigen::MatrixXd V = svd.matrixV();

		Eigen::MatrixXd S = Eigen::MatrixXd::Identity(m, m);

		double det = sigma.determinant();
		double detU = U.determinant();
		double detV = V.determinant();

		Eigen::FullPivLU<Eigen::MatrixXd> luDecomp(sigma);

		int rankSigma = luDecomp.rank();

		std::cout << "rank = " << rankSigma << std::endl;

		if (det < 0 || (rankSigma == m - 1 && detU*detV < 0))
		{
			S(m - 1, m - 1) = -1;
		}

		//Bootstrap

		Eigen::Matrix3d R = U * S*V.transpose();

		Eigen::Vector3d t = mean_Y - R * mean_X;

		poseR = R;
		poseT = t;


	}// end of function

}//end of namespace
