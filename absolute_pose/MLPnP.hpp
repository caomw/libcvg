
#ifndef __CDLT_HPP__
#define __CDLT_HPP__


#pragma warning(disable:4996)

#include <Eigen/Eigen>
#include <Eigen/Core>

namespace cvg
{


	class CMLPnP
	{

	public:
		CMLPnP();
		virtual ~CMLPnP();


		/***************************************/
		/* mR表示旋转矩阵.
		/* mT 表示平移向量.
		/* matXYZ 表示三维点.
		/* matxy 表示二维点.
		/* tol 表示误差容忍度.
		/**************************************/
		int ComputeRT(
			Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
			Eigen::Matrix3Xd &MatXYZ, Eigen::Matrix2Xd &Matxy, double tol = 1e-8);

	private:


		/******************************/
		/* 向量转换为旋转矩阵.
		/* omega表示输出向量.
		/* RotM表示输出举证.
		/******************************/

		int  RodrigueToRotM(Eigen::Vector3d &omega, Eigen::Matrix3d &RotM);


		/******************************/
		/* 矩阵转换为向量.
		/* RotM表示输入向量.
		/* omega表示输出向量.
		/******************************/

		int RodrigueToVect(Eigen::Matrix3d &RotM, Eigen::Vector3d &omega);


		/******************************/
		/* 计算Jacobian矩阵.
		/* jacM 表示输出的Jacobian矩阵, 结果为2行x6列的矩阵.
		/*
		/******************************/

		int JacobiansRodrigues(
			double X1, double Y1, double Z1,
			double r1, double r2, double r3,
			double s1, double s2, double s3,
			double t1, double t2, double t3,
			double w1, double w2, double w3, Eigen::MatrixXd &jacM);


		/******************************/
		/* 计算Jacobian矩阵.
		/* x表示6个元素的列向量.
		/* r表示3行x列的矩阵.
		/* s表示3行x列的矩阵.
		/* points3D 表示3行x列的矩阵(存储世界坐标系下的三维点).
		/* err 表示输出的误差.
		/* J 表示输出的Jacobian矩阵.
		/******************************/

		int ResidualsAndJacobian(
			Eigen::VectorXd &x, Eigen::Matrix3Xd &r, Eigen::Matrix3Xd &s,
			Eigen::Matrix3Xd &points3D, Eigen::VectorXd &err, Eigen::MatrixXd &J);


		struct Statistics
		{
			double resV;
			Eigen::VectorXd r;
			Eigen::Matrix<double, Eigen::Dynamic, 6, Eigen::ColMajor> J;
			Eigen::Matrix<double, 6, 6, Eigen::ColMajor> Qxx;
			double s0;
			Eigen::MatrixXd Qldld;
			Eigen::Matrix<double, 6, 1, Eigen::ColMajor> param;

		};


		struct OptimFlags
		{

		public:
			double epsP;
			double epsF;
			int maxit;
			double tau;

		public:

			OptimFlags()
			{
				epsP = 1e-6;
				epsF = 1e-6;
				maxit = 5;
				tau = 1e-4;
			}

		};


		/******************************************************/
		/* 计算Jacobian矩阵.
		/* Tinit, 输入参数, 表示3行x4列矩阵.
		/* points3D, 输入参数, 表示世界坐标系下的三维点.
		/* rnull, 输入参数,
		/* snull, 输入参数.
		/* P, 输入参数, 表示2n x 2n 矩阵(n表示输入点数).
		/* optimFlags, 输入参数.
		/* T, 输出参数, 表示3x4矩阵.
		/* statistics, 输出参数,
		/*******************************************************/


		int OptimMLPnP_GN(
			Eigen::Matrix<double, 3, 4, Eigen::ColMajor> &Tinit,
			Eigen::Matrix3Xd &points3D, Eigen::Matrix3Xd &rnull,
			Eigen::Matrix3Xd &snull, Eigen::MatrixXd &P,
			const OptimFlags &optimFlags,
			Eigen::Matrix<double, 3, 4, Eigen::ColMajor> &Tout, Statistics &statistics);



		/******************************************************/
		/* MLPnP solver.
		/* points3D, 输入参数, 表示世界坐标系下的三维点.
		/* v, 输入参数, 图像点(齐次坐标表示)
		/* T, 输出参数, 表示3x4矩阵.
		/* statistics, 输出参数,
		/*******************************************************/
		/*Eigen::Matrix<double, 3, 4, Eigen::ColMajor> &T*/

		int MLPnPSolver_Without_Cov(
			Eigen::Matrix3Xd &points3D, Eigen::Matrix3Xd &v,
			Eigen::MatrixXd &T, Statistics &statistics);


		/******************************************************/
		/* MLPnP solver with COV.
		/* points3D, 输入参数, 表示世界坐标系下的三维点.
		/* v, 输入参数, 图像点(齐次坐标表示)
		/* cov, 输入参数,表示协方差矩阵.
		/* T, 输出参数, 表示3x4矩阵.
		/* statistics, 输出参数,
		/*******************************************************/

		int MLPnPSolver_COV(
			Eigen::Matrix3Xd &points3D, Eigen::Matrix3Xd &v,
			Eigen::MatrixXd &cov, Eigen::MatrixXd &T, Statistics &statistics);


	};

}//namespace cvg;


#endif
