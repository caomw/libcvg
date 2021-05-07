

#pragma warning(disable:4996)


#include "LHM.hpp"

#include <iostream>
#include <vector>


namespace cvg
{

	CLHM::CLHM()
	{

	}


	CLHM::~CLHM()
	{

	}



	int CLHM::ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
		Eigen::Matrix3Xd &matXYZ, Eigen::Matrix2Xd &matXY, double tol /*= 1e-8*/)
	{
		if (matXYZ.cols() == 0 || matXYZ.rows() == 0)
		{
			std::cout << "matXYZ is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		if (matXY.cols() == 0 || matXY.rows() == 0)
		{
			std::cout << "matXY is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		//

		double EPSILON = 1e-8;

		std::size_t n = matXYZ.cols();

		Eigen::Vector3d  pbar = matXYZ.rowwise().sum() / (n * 1.0);

		matXYZ = matXYZ.colwise() - pbar;

		Eigen::Matrix3Xd Q = Eigen::Matrix3Xd::Ones(3, n);

		Q.topRows<2>() = matXY;

		//compute projection matrices;

		std::vector<Eigen::Matrix3d> F;

		Eigen::Vector3d V;

		for (std::size_t i = 0; i < n; i++)
		{
			V = Q.col(i) / Q(2, i);
			Eigen::Matrix3d temp = (V* V.transpose()) / (V.transpose()* V);

			F.push_back(temp);

		}

		Eigen::Matrix3d temp2 = Eigen::Matrix3d::Zero();

		for (std::size_t i = 0; i < n; i++)
		{
			temp2 += F[i];

		}

		temp2 = temp2 / (1.0 * n);

		Eigen::Matrix3d tFactor = (Eigen::Matrix3d::Identity() - temp2).inverse() / (1.0 * n);

		//
		//------------------------------

		Eigen::Matrix3d Ri;
		Eigen::Vector3d ti;
		Eigen::Matrix3Xd Qi, Qi2;

		double old_err = 0;
		std::size_t it = 1;

		absKernel(Ri, ti, Qi2, old_err, matXYZ, Q, F, tFactor);

		double new_err = 0.0;

		Qi = Qi2;

		absKernel(Ri, ti, Qi2, new_err, matXYZ, Qi, F, tFactor);
		it++;

		while (std::abs(old_err - new_err) / old_err > tol && new_err > EPSILON)
		{

			old_err = new_err;

			Qi = Qi2;
			absKernel(Ri, ti, Qi2, new_err, matXYZ, Qi, F, tFactor);

			it++;
			if (it > 20)
			{
				break;
			}


		}

		mR = Ri;
		mT = ti;

		mT = mT - Ri * pbar;



		return 1;

	}

	int CLHM::qmatQ(const Eigen::Vector4d &vecQ, Eigen::Matrix4d &matQ)
	{
		double w = vecQ(0);
		double x = vecQ(1);
		double y = vecQ(2);
		double z = vecQ(3);

		matQ(0, 0) = w;
		matQ(0, 1) = -x;
		matQ(0, 2) = -y;
		matQ(0, 3) = -z;

		matQ(1, 0) = x;
		matQ(1, 1) = w;
		matQ(1, 2) = -z;
		matQ(1, 3) = y;

		matQ(2, 0) = y;
		matQ(2, 1) = z;
		matQ(2, 2) = w;
		matQ(2, 3) = -x;

		matQ(3, 0) = z;
		matQ(3, 1) = -y;
		matQ(3, 2) = x;
		matQ(3, 3) = w;

		return 1;


	}

	int CLHM::qmatW(const Eigen::Vector4d &vecQ, Eigen::Matrix4d &matW)
	{
		double w = vecQ(0);
		double x = vecQ(1);
		double y = vecQ(2);
		double z = vecQ(3);

		matW(0, 0) = w;
		matW(0, 1) = -x;
		matW(0, 2) = -y;
		matW(0, 3) = -z;

		matW(1, 0) = x;
		matW(1, 1) = w;
		matW(1, 2) = z;
		matW(1, 3) = -y;

		matW(2, 0) = y;
		matW(2, 1) = -z;
		matW(2, 2) = w;
		matW(2, 3) = x;

		matW(3, 0) = z;
		matW(3, 1) = y;
		matW(3, 2) = -x;
		matW(3, 3) = w;

		return 1;

	}

	int CLHM::quat3Mat(const Eigen::Vector4d &vecQ, Eigen::Matrix3d &matR)
	{
		double a = vecQ(0);
		double b = vecQ(1);
		double c = vecQ(2);
		double d = vecQ(3);

		matR(0, 0) = a * a + b * b - c * c - d * d;
		matR(0, 1) = 2 * (b*c - a * d);
		matR(0, 2) = 2 * (b*d + a * c);

		matR(1, 0) = 2 * (b*c + a * d);
		matR(1, 1) = a * a + c * c - b * b - d * d;
		matR(1, 2) = 2 * (c*d - a * b);

		matR(2, 0) = 2 * (b*d - a * c);
		matR(2, 1) = 2 * (c*d + a * b);
		matR(2, 2) = a * a + d * d - b * b - c * c;

		return 1;

	}

	int CLHM::xForm(
		const Eigen::Matrix3Xd &matP,
		const Eigen::Matrix3d &matR,
		const Eigen::Vector3d &vecT,
		Eigen::Matrix3Xd &matQ)
	{
		std::size_t n = matP.cols();

		if (n == 0)
		{
			std::cout << "matP is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;
			return -1;
		}

		matQ = Eigen::Matrix3Xd::Zero(3, n);

		for (std::size_t i = 0; i < n; i++)
		{
			matQ.col(i) = matR * matP.col(i) + vecT;

		}

		return 1;

	}

	int CLHM::xFormProj(const Eigen::Matrix3Xd &matP, const Eigen::Matrix3d &matR,
		const Eigen::Vector3d &vecT, Eigen::Matrix2Xd &matQp)
	{

		std::size_t n = matP.cols();
		if (n == 0)
		{
			std::cout << "matP is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}

		Eigen::Matrix3Xd matQ = Eigen::Matrix3Xd::Zero(3, n);

		for (std::size_t i = 0; i < n; i++)
		{
			matQ.col(i) = matR * matP.col(i) + vecT;

		}

		matQp = Eigen::Matrix2Xd::Zero(2, n);

		matQp = matQ.topRows<2>();
		matQp.row(0) = matQp.row(0).cwiseQuotient(matQ.row(3));
		matQp.row(1) = matQp.row(1).cwiseQuotient(matQ.row(3));

		return 1;

	}

	int CLHM::estimate_t(Eigen::Vector3d &t, const Eigen::Matrix3d &R, const Eigen::Matrix3d &G,
		const std::vector<Eigen::Matrix3d> &F, const Eigen::Matrix3Xd &P)
	{
		std::size_t n = P.cols();

		if (n == 0)
		{
			std::cout << "P is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;
			return -1;
		}

		Eigen::Vector3d sum_ = Eigen::Vector3d::Zero();

		for (std::size_t i = 0; i < n; i++)
		{
			sum_ = sum_ + F[i] * R * P.col(i);

		}

		//计算返回值.
		t = G * sum_;


		return 1;

	}

	int CLHM::absKernel(Eigen::Matrix3d &R, Eigen::Vector3d &t, Eigen::Matrix3Xd &Qout,
		double &err2, const Eigen::Matrix3Xd &P, const Eigen::Matrix3Xd &Q,
		const std::vector<Eigen::Matrix3d> &F, const Eigen::Matrix3d &G)
	{
		std::size_t n = P.cols();

		if (n == 0)
		{
			std::cout << "P is empty !" << " at line: " <<
				__LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		Eigen::Matrix3Xd matQ = Q;

		for (std::size_t i = 0; i < n; i++)
		{
			matQ.col(i) = F[i] * matQ.col(i);

		}

		Eigen::Vector3d pbar = P.rowwise().sum() / (1.0 * n);


		Eigen::Matrix3Xd matP = P;

		Eigen::Vector3d qbar = matQ.rowwise().sum() / (1.0 * n);


		matP = matP.colwise() - pbar;

		matQ = matQ.colwise() - qbar;


		Eigen::Matrix3d M = Eigen::Matrix3d::Zero();

		for (std::size_t i = 0; i < n; i++)
		{
			M += (matP.col(i) * matQ.col(i).transpose());
		}

		Eigen::JacobiSVD<Eigen::Matrix3d> svd(M, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::Matrix3d U = svd.matrixU();
		Eigen::Matrix3d V = svd.matrixV();

		R = V * U.transpose();

		if (R.determinant() > 0)
		{

			estimate_t(t, R, G, F, matP);

			if (t(2) < 0)
			{
				V.col(2) = -V.col(2);
				R = -(V * U.transpose());

				estimate_t(t, R, G, F, matP);
			}


		}
		else
		{

			V.col(2) = -V.col(2);

			R = V * U.transpose();

			estimate_t(t, R, G, F, matP);

			if (t(2) < 0)
			{

				R = -V * U.transpose();

				estimate_t(t, R, G, F, matP);

			}


		}


		xForm(matP, R, t, Qout);

		err2 = 0;

		for (std::size_t i = 0; i < n; i++)
		{
			Eigen::Vector3d vec = (Eigen::Matrix3d::Identity() - F[i])*Qout.col(i);

			err2 += (vec.dot(vec));

		}

		return 1;

	}

}
