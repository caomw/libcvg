
#ifndef __CDLT_HPP__
#define __CDLT_HPP__


#pragma warning(disable:4996)

#include <Eigen/Eigen>
#include <Eigen/Core>


namespace cvg
{


	class CDLT
	{
	public:
		CDLT();
		virtual ~CDLT();


		/***************************************/
		/* mR表示旋转矩阵.
		/* mT 表示平移向量.
		/* matXYZ 表示三维点.
		/* matxy 表示二维点.
		/* tol 表示误差容忍度.
		/**************************************/
		int ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
			Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, double tol = 1e-8);

	private:
		int CalcPose(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
			Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix3Xd &Matxy);


		/*******************************/
		/* 符号函数.

		/*******************************/

		int sign(double x);


		int EigenSolver(Eigen::MatrixXd &mat, Eigen::MatrixXd &matV);


	};

}//namespace cvg;





#endif
