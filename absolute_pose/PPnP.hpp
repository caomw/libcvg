#pragma once

#pragma warning(disable:4996)

#include <string>


#include <Eigen/Core>
#include <Eigen/Dense>



class CPPnP
{
public:
	CPPnP();
	virtual ~CPPnP();


public:

	/*
		/ *
		* pName 表示二维点文件名称.
		* sName 表示三维点文件名称.
		* /
		int LoadData(const std::string pName, const std::string sName);


		/ *
		* mR 表示旋转矩阵.
		* mT 表示平移向量.
		* tol 表示误差容忍度.
		* /
		int ComputeRT(Eigen::Matrix3d & mR, Eigen::Vector3d &mT, const double tol = 1e-8);
		*/


		/*
		   * mR 表示旋转矩阵.
		   * mT 表示平移向量.
		   * tol 表示误差容忍度.
		   * mP 表示输入二维点(采用齐次坐标表示)
		   * mS 表示世界坐标系中的三维点.
		   */
	int ComputeRT(
		Eigen::Matrix3d & mR, Eigen::Vector3d &mT,
		const Eigen::MatrixXd &tmS, const Eigen::MatrixXd &tmP,
		const double tol = 1e-8);

private:

	/*
	  *mP 表示输入二维点(采用齐次坐标表示)
	  * mS 表示世界坐标系中的三维点.
	  */
	Eigen::MatrixXd mP;

	Eigen::MatrixXd mS;


};
