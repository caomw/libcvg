

#pragma warning(disable:4996)

#include "P3P.hpp"


#include <cmath>
#include <algorithm>
#include <iostream>



CP3P::CP3P()
{

}


CP3P::~CP3P()
{

}

int CP3P::ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
	Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, unsigned int index/* =0*/)
{
	if (MatXYZ.cols() < 3)
	{
		std::cout << "MatXYZ is empty !" << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;

		return -1;
	}


	if (Matxy.cols() < 3)
	{
		std::cout << "Matxy is empty !" << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;

		return -1;
	}

	if (index >= 4 || index < 0)
	{
		std::cout << "Array out of index !" << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;

		return -1;
	}


	std::vector<Eigen::Matrix3d> matR;
	std::vector<Eigen::Vector3d> matT;

	if (KComputePoses(matR, matT, MatXYZ, Matxy))
	{
		mR = matR[index];
		mT = matT[index];
	}


	return 1;

}




int CP3P::ComputePoses(Eigen::Matrix<double, 3, 16, Eigen::ColMajor> &poses,
	Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy)
{

	if (MatXYZ.cols() < 3 || Matxy.cols() < 3)
	{
		std::cout << "The input matrix are incorrect !" << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;

		return -1;
	}

	//取世界坐标系内的前3个点.
	Eigen::Matrix3d wXYZ;

	if (MatXYZ.cols() > 3)
	{
		wXYZ = MatXYZ.leftCols<3>();
	}
	else
	{
		wXYZ = MatXYZ;
	}


	Eigen::Vector3d P1 = wXYZ.col(0);
	Eigen::Vector3d P2 = wXYZ.col(1);
	Eigen::Vector3d P3 = wXYZ.col(2);

	//判断世界坐标系内的点是否共面.
	Eigen::Vector3d vector1 = P2 - P1;
	Eigen::Vector3d vector2 = P3 - P1;

	if ((vector1.cross(vector2)).norm() == 2)
	{
		std::cout << "The inputed points are coplanar !" << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;

		return -1;
	}

	//取图像坐标系下的前3个点.
	Eigen::Matrix3d matTemp = Eigen::Matrix3d::Ones();

	if (Matxy.cols() > 3)
	{
		matTemp.topRows<2>() = Matxy.leftCols<3>();
	}
	else
	{
		matTemp.topRows<2>() = Matxy;
	}


	Eigen::Vector3d f1 = matTemp.col(0);
	Eigen::Vector3d f2 = matTemp.col(1);
	Eigen::Vector3d f3 = matTemp.col(2);

	Eigen::Vector3d e1 = f1;
	Eigen::Vector3d e3 = f1.cross(f2);

	double norm_e3 = e3.norm();

	Eigen::Vector3d temp;
	temp.fill(norm_e3);

	e3 = e3.cwiseQuotient(temp);

	Eigen::Vector3d e2 = e3.cross(e1);

	Eigen::Matrix3d matT;

	matT.row(0) = e1.transpose();
	matT.row(1) = e2.transpose();
	matT.row(2) = e3.transpose();

	f3 = matT * f3;

	if (f3(2) > 0)
	{
		Eigen::Vector3d f1 = matTemp.col(1);
		Eigen::Vector3d f2 = matTemp.col(0);
		Eigen::Vector3d f3 = matTemp.col(2);

		e1 = f1;
		e3 = f1.cross(f2);

		norm_e3 = e3.norm();
		temp.fill(norm_e3);
		e3.cwiseQuotient(temp);

		e2 = e3.cross(e1);

		matT.row(0) = e1.transpose();
		matT.row(1) = e2.transpose();
		matT.row(2) = e3.transpose();

		f3 = matT * f3;

		P1 = wXYZ.col(1);
		P2 = wXYZ.col(0);
		P3 = wXYZ.col(2);

	}


	Eigen::Vector3d n1 = P2 - P1;

	double norm_n1 = n1.norm();

	Eigen::Vector3d temp1;

	temp1.fill(norm_n1);
	n1 = n1.cwiseQuotient(temp1);

	Eigen::Vector3d n3 = n1.cross((P3 - P1));


	double norm_n3 = n3.norm();

	Eigen::Vector3d temp2;
	temp2.fill(norm_n3);

	n3 = n3.cwiseQuotient(temp2);

	Eigen::Vector3d n2 = n3.cross(n1);

	Eigen::Matrix3d matN;
	matN.row(0) = n1.transpose();
	matN.row(1) = n2.transpose();
	matN.row(2) = n3.transpose();

	//已知参数提取.
	P3 = matN * (P3 - P1);

	double d_12 = (P2 - P1).norm();
	double f_1 = f3(0) / f3(2);
	double f_2 = f3(1) / f3(2);
	double p_1 = P3(0);
	double p_2 = P3(1);

	double cos_beta = f1.transpose()*f2;

	double b = 1.0 / (1 - cos_beta * cos_beta) - 1;

	if (cos_beta < 0)
	{
		b = -std::sqrt(b);
	}
	else
	{
		b = std::sqrt(b);
	}


	//定义一些变量避免重复计算.

	double f_1_pw2 = f_1 * f_1;
	double f_2_pw2 = f_2 * f_2;
	double p_1_pw2 = p_1 * p_1;
	double p_1_pw3 = p_1_pw2 * p_1;
	double p_1_pw4 = p_1_pw3 * p_1;
	double p_2_pw2 = p_2 * p_2;
	double p_2_pw3 = p_2_pw2 * p_2;
	double p_2_pw4 = p_2_pw3 * p_2;
	double d_12_pw2 = d_12 * d_12;
	double b_pw2 = b * b;




	//Computation of factors of 4th degree polynomial

	double factor_4 = -f_2_pw2 * p_2_pw4
		- p_2_pw4 * f_1_pw2 - p_2_pw4;

	double factor_3 = 2 * p_2_pw3*d_12*b
		+ 2 * f_2_pw2*p_2_pw3*d_12*b - 2 * f_2*p_2_pw3*f_1*d_12;

	double factor_2 = -f_2_pw2 * p_2_pw2*p_1_pw2
		- f_2_pw2 * p_2_pw2*d_12_pw2*b_pw2
		- f_2_pw2 * p_2_pw2*d_12_pw2
		+ f_2_pw2 * p_2_pw4
		+ p_2_pw4 * f_1_pw2
		+ 2 * p_1*p_2_pw2*d_12
		+ 2 * f_1*f_2*p_1*p_2_pw2*d_12*b
		- p_2_pw2 * p_1_pw2*f_1_pw2
		+ 2 * p_1*p_2_pw2*f_2_pw2*d_12
		- p_2_pw2 * d_12_pw2*b_pw2
		- 2 * p_1_pw2*p_2_pw2;

	double factor_1 = 2 * p_1_pw2*p_2*d_12*b
		+ 2 * f_2*p_2_pw3*f_1*d_12
		- 2 * f_2_pw2*p_2_pw3*d_12*b
		- 2 * p_1*p_2*d_12_pw2*b;

	double factor_0 = -2 * f_2*p_2_pw2*f_1*p_1*d_12*b
		+ f_2_pw2 * p_2_pw2*d_12_pw2
		+ 2 * p_1_pw3*d_12
		- p_1_pw2 * d_12_pw2
		+ f_2_pw2 * p_2_pw2*p_1_pw2
		- p_1_pw4
		- 2 * f_2_pw2*p_2_pw2*p_1*d_12
		+ p_2_pw2 * f_1_pw2*p_1_pw2
		+ f_2_pw2 * p_2_pw2*d_12_pw2*b_pw2;


	Eigen::VectorXd factors = Eigen::VectorXd::Zero(5);

	factors(0) = factor_4;
	factors(1) = factor_3;
	factors(2) = factor_2;
	factors(3) = factor_1;
	factors(4) = factor_0;

	//Computation of roots

	Eigen::Vector4d x;


	SolveQuartic(factors, x);



	//Backsubstitution of each solution


	for (int i = 0; i < 4; i++)
	{


		double cot_alpha = (-f_1 * p_1 / f_2 - x(i)*p_2 + d_12 * b) / (-f_1 * x(i)*p_2 / f_2 + p_1 - d_12);

		double cos_theta = x(i);
		double sin_theta = std::sqrt(1.0 - x(i)*x(i));
		double sin_alpha = std::sqrt(1.0 / (cot_alpha * cot_alpha + 1));
		double cos_alpha = std::sqrt(1.0 - sin_alpha * sin_alpha);



		if (cot_alpha < 0)
		{
			cos_alpha = -cos_alpha;
		}

		Eigen::Vector3d vecC;
		vecC(0) = d_12 * cos_alpha*(sin_alpha*b + cos_alpha);
		vecC(1) = cos_theta * d_12*sin_alpha*(sin_alpha*b + cos_alpha);
		vecC(2) = sin_theta * d_12*sin_alpha*(sin_alpha*b + cos_alpha);

		vecC = P1 + matN.transpose()* vecC;


		Eigen::Matrix3d matR;
		matR(0, 0) = -cos_alpha;
		matR(0, 1) = -sin_alpha * cos_theta;
		matR(0, 2) = -sin_alpha * sin_theta;

		matR(1, 0) = sin_alpha;
		matR(1, 1) = -cos_alpha * cos_theta;
		matR(1, 2) = -cos_alpha * sin_theta;

		matR(2, 0) = 0;
		matR(2, 1) = -sin_theta;
		matR(2, 2) = cos_theta;


		matR = matN.transpose() * matR.transpose() * matT;

		poses.col(i * 4) = vecC;

		//poses.middleCols<3>(i * 4 + 1) = matR;

		poses.col(i * 4 + 1) = matR.col(0);

		poses.col(i * 4 + 2) = matR.col(1);

		poses.col(i * 4 + 3) = matR.col(2);

	}

	return -1;
}

int CP3P::KComputePoses(
	std::vector<Eigen::Matrix3d> &mR,
	std::vector<Eigen::Vector3d>&mT,
	Eigen::Matrix3Xd &MatXYZ,
	Eigen::Matrix2Xd &Matxy)
{

	if (MatXYZ.cols() < 3 || Matxy.cols() < 3)
	{
		std::cout << "The input matrices are incorrect !"
			<< " at line: " << __LINE__ << " in file: " << __FILE__ << std::endl;
		return -1;

	}

	Eigen::Matrix3d wXYZ;


	if (MatXYZ.cols() > 3)
	{
		wXYZ = MatXYZ.leftCols<3>();
	}
	else
	{
		wXYZ = MatXYZ;
	}

	Eigen::Matrix3d matImgPt = Eigen::Matrix3d::Ones();

	if (Matxy.cols() > 3)
	{
		matImgPt.topRows<2>() = Matxy.leftCols<3>();
	}
	else
	{
		matImgPt.topRows<2>() = Matxy;
	}

	Eigen::Matrix3d matTemp = matImgPt.cwiseProduct(matImgPt);
	Eigen::RowVector3d vecTemp = matTemp.colwise().sum();


	vecTemp = vecTemp.cwiseSqrt();

	matTemp.row(0) = vecTemp;
	matTemp.row(1) = vecTemp;
	matTemp.row(2) = vecTemp;

	matImgPt = matImgPt.cwiseQuotient(matTemp);

	////////////////////

	//Eigen::Matrix<double, 3, 16, Eigen::ColMajor> poses;

	Eigen::Vector3d P1 = wXYZ.col(0);
	Eigen::Vector3d P2 = wXYZ.col(1);
	Eigen::Vector3d P3 = wXYZ.col(2);

	//判断世界坐标系内的点是否共面.
	Eigen::Vector3d vector1 = P2 - P1;
	Eigen::Vector3d vector2 = P3 - P1;

	if ((vector1.cross(vector2)).norm() == 2)
	{
		std::cout << "The inputed points are coplanar !" << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;

		return -1;
	}


	Eigen::Vector3d f1 = matImgPt.col(0);
	Eigen::Vector3d f2 = matImgPt.col(1);
	Eigen::Vector3d f3 = matImgPt.col(2);

	Eigen::Vector3d e1 = f1;
	Eigen::Vector3d e3 = f1.cross(f2);

	double norm_e3 = e3.norm();

	Eigen::Vector3d temp;
	temp.fill(norm_e3);

	e3 = e3.cwiseQuotient(temp);

	Eigen::Vector3d e2 = e3.cross(e1);

	Eigen::Matrix3d matT;

	matT.row(0) = e1.transpose();
	matT.row(1) = e2.transpose();
	matT.row(2) = e3.transpose();

	f3 = matT * f3;

	if (f3(2) > 0)
	{
		Eigen::Vector3d f1 = matImgPt.col(1);
		Eigen::Vector3d f2 = matImgPt.col(0);
		Eigen::Vector3d f3 = matImgPt.col(2);

		e1 = f1;
		e3 = f1.cross(f2);

		norm_e3 = e3.norm();

		Eigen::Vector3d temp;
		temp.fill(norm_e3);
		e3 = e3.cwiseQuotient(temp);

		e2 = e3.cross(e1);

		matT.row(0) = e1.transpose();
		matT.row(1) = e2.transpose();
		matT.row(2) = e3.transpose();

		f3 = matT * f3;

		P1 = wXYZ.col(1);
		P2 = wXYZ.col(0);
		P3 = wXYZ.col(2);

	}


	Eigen::Vector3d n1 = P2 - P1;

	double norm_n1 = n1.norm();

	Eigen::Vector3d temp1;

	temp1.fill(norm_n1);
	n1 = n1.cwiseQuotient(temp1);

	Eigen::Vector3d n3 = n1.cross(P3 - P1);


	double norm_n3 = n3.norm();

	Eigen::Vector3d temp2;
	temp2.fill(norm_n3);

	n3 = n3.cwiseQuotient(temp2);

	Eigen::Vector3d n2 = n3.cross(n1);

	Eigen::Matrix3d matN;
	matN.row(0) = n1.transpose();
	matN.row(1) = n2.transpose();
	matN.row(2) = n3.transpose();

	//已知参数提取.
	P3 = matN * (P3 - P1);

	double d_12 = (P2 - P1).norm();
	double f_1 = f3(0) / f3(2);
	double f_2 = f3(1) / f3(2);
	double p_1 = P3(0);
	double p_2 = P3(1);

	double cos_beta = f1.transpose()*f2;

	double b = 1.0 / (1.0 - cos_beta * cos_beta) - 1;

	if (cos_beta < 0)
	{
		b = -std::sqrt(b);
	}
	else
	{
		b = std::sqrt(b);
	}


	//定义一些变量避免重复计算.

	double f_1_pw2 = f_1 * f_1;
	double f_2_pw2 = f_2 * f_2;
	double p_1_pw2 = p_1 * p_1;
	double p_1_pw3 = p_1_pw2 * p_1;
	double p_1_pw4 = p_1_pw3 * p_1;
	double p_2_pw2 = p_2 * p_2;
	double p_2_pw3 = p_2_pw2 * p_2;
	double p_2_pw4 = p_2_pw3 * p_2;
	double d_12_pw2 = d_12 * d_12;
	double b_pw2 = b * b;




	//Computation of factors of 4th degree polynomial

	double factor_4 = -f_2_pw2 * p_2_pw4
		- p_2_pw4 * f_1_pw2 - p_2_pw4;

	double factor_3 = 2 * p_2_pw3*d_12*b
		+ 2 * f_2_pw2*p_2_pw3*d_12*b - 2 * f_2*p_2_pw3*f_1*d_12;

	double factor_2 = -f_2_pw2 * p_2_pw2*p_1_pw2
		- f_2_pw2 * p_2_pw2*d_12_pw2*b_pw2
		- f_2_pw2 * p_2_pw2*d_12_pw2
		+ f_2_pw2 * p_2_pw4
		+ p_2_pw4 * f_1_pw2
		+ 2 * p_1*p_2_pw2*d_12
		+ 2 * f_1*f_2*p_1*p_2_pw2*d_12*b
		- p_2_pw2 * p_1_pw2*f_1_pw2
		+ 2 * p_1*p_2_pw2*f_2_pw2*d_12
		- p_2_pw2 * d_12_pw2*b_pw2
		- 2 * p_1_pw2*p_2_pw2;

	double factor_1 = 2 * p_1_pw2*p_2*d_12*b
		+ 2 * f_2*p_2_pw3*f_1*d_12
		- 2 * f_2_pw2*p_2_pw3*d_12*b
		- 2 * p_1*p_2*d_12_pw2*b;

	double factor_0 = -2 * f_2*p_2_pw2*f_1*p_1*d_12*b
		+ f_2_pw2 * p_2_pw2*d_12_pw2
		+ 2 * p_1_pw3*d_12
		- p_1_pw2 * d_12_pw2
		+ f_2_pw2 * p_2_pw2*p_1_pw2
		- p_1_pw4
		- 2 * f_2_pw2*p_2_pw2*p_1*d_12
		+ p_2_pw2 * f_1_pw2*p_1_pw2
		+ f_2_pw2 * p_2_pw2*d_12_pw2*b_pw2;


	Eigen::VectorXd factors = Eigen::VectorXd::Zero(5);

	factors(0) = factor_4;
	factors(1) = factor_3;
	factors(2) = factor_2;
	factors(3) = factor_1;
	factors(4) = factor_0;

	//Computation of roots

	Eigen::Vector4d x;


	SolveQuartic(factors, x);



	//Backsubstitution of each solution

	Eigen::Matrix3d matR;
	Eigen::Vector3d vecT;

	for (int i = 0; i < 4; i++)
	{


		double cot_alpha = (-f_1 * p_1 / f_2 - x(i)*p_2 + d_12 * b) / (-f_1 * x(i)*p_2 / f_2 + p_1 - d_12);

		double cos_theta = x(i);
		double sin_theta = std::sqrt(1.0 - x(i)*x(i));
		double sin_alpha = std::sqrt(1.0 / (cot_alpha * cot_alpha + 1));
		double cos_alpha = std::sqrt(1.0 - sin_alpha * sin_alpha);



		if (cot_alpha < 0)
		{
			cos_alpha = -cos_alpha;
		}

		Eigen::Vector3d vecC;
		vecC(0) = d_12 * cos_alpha*(sin_alpha*b + cos_alpha);
		vecC(1) = cos_theta * d_12*sin_alpha*(sin_alpha*b + cos_alpha);
		vecC(2) = sin_theta * d_12*sin_alpha*(sin_alpha*b + cos_alpha);

		vecC = P1 + matN.transpose()* vecC;

		matR(0, 0) = -cos_alpha;
		matR(0, 1) = -sin_alpha * cos_theta;
		matR(0, 2) = -sin_alpha * sin_theta;

		matR(1, 0) = sin_alpha;
		matR(1, 1) = -cos_alpha * cos_theta;
		matR(1, 2) = -cos_alpha * sin_theta;

		matR(2, 0) = 0;
		matR(2, 1) = -sin_theta;
		matR(2, 2) = cos_theta;

		matR = matN.transpose() * matR.transpose() * matT;


		vecT = -matR.transpose() * vecC;

		Eigen::Matrix3d matProj = matR.transpose() *wXYZ;

		matProj = matProj.colwise() + vecT;

		if (matProj.row(2).minCoeff() > 0)
		{
			mR.push_back(matR.transpose());
			mT.push_back(vecT);
		}


	}

	//////////////////////////////


	return 1;

}

int CP3P::SolveQuartic(Eigen::VectorXd &factors, Eigen::Vector4d &roots)
{

	if (factors.size() != 5)
	{
		std::cout << "factors.size()!=5 ." << " at line: "
			<< __LINE__ << " in file: " << __FILE__ << std::endl;
		return -1;
	}


	double A = factors(0);
	double B = factors(1);
	double C = factors(2);
	double D = factors(3);
	double E = factors(4);

	double A_pw2 = A * A;
	double B_pw2 = B * B;
	double A_pw3 = A_pw2 * A;
	double B_pw3 = B_pw2 * B;
	double A_pw4 = A_pw3 * A;
	double B_pw4 = B_pw3 * B;

	double alpha = -3.0 * B_pw2 / (8.0 * A_pw2) + C / A;
	double beta = B_pw3 / (8.0 * A_pw3) - B * C / (2.0 * A_pw2) + D / A;
	double gamma = -3.0 * B_pw4 / (256.0 * A_pw4) + B_pw2 * C / (16.0 * A_pw3) - B * D / (4.0 * A_pw2) + E / A;

	double alpha_pw2 = alpha * alpha;
	double alpha_pw3 = alpha_pw2 * alpha;

	std::complex<double> P(-alpha_pw2 / 12.0 - gamma, 0);
	std::complex<double> Q(-alpha_pw3 / 108.0 + alpha * gamma / 3.0 - beta * beta / 8.0, 0);

	std::complex<double> R = -Q / 2.0 + std::sqrt(Q*Q / 4.0 + P * P*P / 27.0);

	std::complex<double> U = std::pow(R, 1.0 / 3.0);

	std::complex<double> y = 0.0;

	if (U.real() == 0)
	{
		y = -5.0 * alpha / 6.0 - std::pow(Q, (1.0 / 3.0));
	}
	else
	{
		y = -5.0*alpha / 6.0 - P / (3.0 * U) + U;
	}


	std::complex<double> w = std::sqrt(alpha + 2.0 *y);

	std::complex<double> temp = 0.0;

	temp = -B / (4.0 * A) + 0.5*(w + std::sqrt(-(3.0 * alpha + 2.0 * y + 2.0 * beta / w)));

	roots(0) = temp.real();



	temp = -B / (4.0 * A) + 0.5*(w - std::sqrt(-(3 * alpha + 2.0 * y + 2.0 * beta / w)));
	roots(1) = temp.real();


	temp = -B / (4.0 * A) + 0.5*(-w + std::sqrt(-(3.0 * alpha + 2.0 * y - 2.0 * beta / w)));
	roots(2) = temp.real();


	temp = -B / (4.0 * A) + 0.5*(-w - std::sqrt(-(3.0 * alpha + 2.0* y - 2.0 * beta / w)));

	roots(3) = temp.real();

	return 1;


}
