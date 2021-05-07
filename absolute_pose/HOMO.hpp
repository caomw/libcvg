
#ifndef __CHOMO_HPP__
#define __CHOMO_HPP__

#pragma warning(disable:4996)

#include <Eigen/Core>


namespace cvg
{


	class CHOMO
	{
	public:
		CHOMO();
		virtual ~CHOMO();

		/************************************************************************/
		/* mR表示旋转矩阵.
		/* mT 表示平移向量.
		/* matXYZ 表示三维点.
		/* matzy 表示二维点.
		/* tol 表示误差容忍度.
		/************************************************************************/
		int ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
			Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, double tol = 1e-8);



		/************************************************************************/
		/* mat表示输入矩阵.
		/* matV 表示矩阵mat的特征向量.
		/************************************************************************/
		int EigenSolver(Eigen::MatrixXd &mat, Eigen::MatrixXd &matV);

		int sign(double x);



	};

}//namespace cvg;




#endif
