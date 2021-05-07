
#pragma warning(disable:4996)


#include "PPnP.hpp"

#include <iostream>

#include <limits>

#include <Eigen/SVD>


CPPnP::CPPnP()
{

}


CPPnP::~CPPnP()
{

}
/*

int CPPnP::LoadData(const std::string pName, const std::string sName)
{

if (pName.empty())
{
std::cout << "pName is empty !" << " at line: "
<< __LINE__ << ", in file: " << __FILE__ << std::endl;

return -1;
}


if (sName.empty())
{
std::cout << "sName is empty !" << " at line: "
<< __LINE__ << ", in file: " << __FILE__ << std::endl;

return -1;
}


//加载数据.




}
*/


/*

int CPPnP::ComputeRT(Eigen::Matrix3d & mR, Eigen::Vector3d &mT, const double tol / *= 1e-8* /)
{
if (tol == 0)
{
std::cout << "tol == 0 !" << " at line: "
<< __LINE__ << ", in file: " << __FILE__ << std::endl;

return -1;
}


std::size_t n = mP.rows();

Eigen::MatrixXd mZ = Eigen::MatrixXd::Zero(n, n);

Eigen::VectorXd mE = Eigen::VectorXd::Ones(n);

Eigen::MatrixXd mA = Eigen::MatrixXd::Identity(n, n) - ((mE * mE.transpose()) / n);

Eigen::VectorXd mII = mE / n;

double err = std::numeric_limits<double>::max();

Eigen::MatrixX3d E_old = 1000 * Eigen::MatrixX3d::Ones(n, 3);

Eigen::Vector3d mC;

//循环开始

while (err>tol)
{

Eigen::Matrix3d mat = mP.transpose()*mZ * mA * mS;

Eigen::JacobiSVD<Eigen::Matrix3d>svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);

Eigen::Matrix3d mU = svd.matrixU();

Eigen::Matrix3d mV = svd.matrixV();

Eigen::Matrix3d mVT = mV.transpose();

Eigen::Matrix3d m = Eigen::Matrix3d::Identity();

m(2, 2) = (mU*mVT).determinant();

mR = mU* m * mVT;

Eigen::MatrixX3d mPR = mP * mR;

mC = (mS - mZ * mPR).transpose() * mII;

Eigen::MatrixX3d mY = mS - mE * mC.transpose();


Eigen::VectorXd Zmindiag = (mPR * mY.transpose()).diagonal().cwiseQuotient(mP.cwiseAbs2().rowwise().sum());

Zmindiag = (Zmindiag.array() < 0).select(0, Zmindiag.array());

mZ = Zmindiag.asDiagonal();

Eigen::MatrixX3d E_new = mY - mZ * mPR;

err = std::sqrt(((E_new - E_old).transpose() * (E_new - E_old)).trace());

E_old = E_new;

}

mT = -mR * mC;


return 1;


}

*/

int CPPnP::ComputeRT(
	Eigen::Matrix3d & mR, Eigen::Vector3d &mT,
	const Eigen::MatrixXd &tmS, const Eigen::MatrixXd &tmP,
	const double tol/*= 1e-8*/)
{


	if (tol == 0)
	{
		std::cout << "tol == 0 !" << " at line: "
			<< __LINE__ << ", in file: " << __FILE__ << std::endl;

		return -1;
	}

	//数据格式转换.
	Eigen::MatrixX3d mS = tmS.transpose();
	Eigen::MatrixX3d mP = Eigen::MatrixX3d::Ones(tmP.cols(), tmP.rows() + 1);

	mP.leftCols<2>() = tmP.transpose();




	std::size_t n = mP.rows();

	Eigen::MatrixXd mZ = Eigen::MatrixXd::Zero(n, n);

	Eigen::VectorXd mE = Eigen::VectorXd::Ones(n);

	Eigen::MatrixXd mA = Eigen::MatrixXd::Identity(n, n) - ((mE * mE.transpose()) / n);

	Eigen::VectorXd mII = mE / n;

	double err = std::numeric_limits<double>::max();

	Eigen::MatrixX3d E_old = 1000 * Eigen::MatrixX3d::Ones(n, 3);

	Eigen::Vector3d mC;

	//循环开始

	while (err > tol)
	{

		Eigen::Matrix3d mat = mP.transpose()*mZ * mA * mS;

		Eigen::JacobiSVD<Eigen::Matrix3d>svd(mat, Eigen::ComputeFullU | Eigen::ComputeFullV);

		Eigen::Matrix3d mU = svd.matrixU();

		Eigen::Matrix3d mV = svd.matrixV();

		Eigen::Matrix3d mVT = mV.transpose();

		Eigen::Matrix3d m = Eigen::Matrix3d::Identity();

		m(2, 2) = (mU*mVT).determinant();

		mR = mU * m * mVT;

		Eigen::MatrixX3d mPR = mP * mR;

		mC = (mS - mZ * mPR).transpose() * mII;

		Eigen::MatrixX3d mY = mS - mE * mC.transpose();


		Eigen::VectorXd Zmindiag = (mPR * mY.transpose()).diagonal().cwiseQuotient(mP.cwiseAbs2().rowwise().sum());

		Zmindiag = (Zmindiag.array() < 0).select(0, Zmindiag.array());

		mZ = Zmindiag.asDiagonal();

		Eigen::MatrixX3d E_new = mY - mZ * mPR;

		err = std::sqrt(((E_new - E_old).transpose() * (E_new - E_old)).trace());

		E_old = E_new;

	}

	mT = -mR * mC;



	return -1;

}
