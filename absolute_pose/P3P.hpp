
#ifndef __CP3P_HPP__
#define __CP3P_HPP__

#pragma warning(disable:4996)


#include <Eigen/Eigen>
#include <Eigen/Core>

#include <vector>


class CP3P
{
public:
	CP3P();
	virtual ~CP3P();


	/***********************************/
	/* mR表示旋转矩阵.
	/* mT 表示平移向量.
	/* matXYZ 表示三维点.
	/* matxy 表示二维点.
	/* index 解的序号.
	/************************************/
	int ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
		Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, unsigned int index = 0);


private:


	/***********************************/
	/* poses表示姿态矩阵,共有4种解.
	/* matXYZ 表示三维点.
	/* matxy 表示二维点.
	/************************************/
	int ComputePoses(Eigen::Matrix<double, 3, 16, Eigen::ColMajor> &poses, Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy);


public:
	/***********************************/
	/* mR表示旋转矩阵集合.
	/* mT表示平移向量集合.
	/* matXYZ 表示三维点.
	/* matxy 表示二维点.
	/************************************/
	int KComputePoses(
		std::vector<Eigen::Matrix3d> &mR,
		std::vector<Eigen::Vector3d>&mT,
		Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy);



	int SolveQuartic(Eigen::VectorXd &factors, Eigen::Vector4d &roots);

private:



};



#endif
