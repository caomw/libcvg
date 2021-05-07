
#ifndef __CDLT_HPP__
#define __CDLT_HPP__

#pragma warning(disable:4996)


#include <Eigen/Eigen>
#include <Eigen/Core>


#include <vector>

namespace cvg
{

	class CLHM
	{
	public:
		CLHM();
		virtual ~CLHM();


		/*********************************/
		/* mR表示旋转矩阵.
		/* mT 表示平移向量.
		/* matXYZ 表示三维点.
		/* matxy 表示二维点.
		/* tol 表示误差容忍度.
		/*********************************/
		int ComputeRT(Eigen::Matrix3d &mR, Eigen::Vector3d &mT,
			Eigen::Matrix3Xd &matXYZ, Eigen::Matrix2Xd &matXY, double tol = 1e-5);

	private:

		/********************************/
		/* 将向量转换为矩阵.
		/*******************************/
		int qmatQ(const Eigen::Vector4d &vecQ, Eigen::Matrix4d &matQ);


		/********************************/
		/* 将向量转换为矩阵.
		/*******************************/
		int qmatW(const Eigen::Vector4d &vecQ, Eigen::Matrix4d &matW);

		/********************************/
		/* 将四元组转换为旋转矩阵.
		/*******************************/
		int quat3Mat(const Eigen::Vector4d &vecQ, Eigen::Matrix3d &matR);


		/********************************/
		/* matP表示世界坐标系下的三维点.
		/* matR表示旋转矩阵.
		/* vecT表示平移向量.
		/* matQ 表示输出结果.
		/*******************************/
		int xForm(const Eigen::Matrix3Xd &matP, const Eigen::Matrix3d &matR,
			const Eigen::Vector3d &vecT, Eigen::Matrix3Xd &matQ);

		/********************************/
		/* matP表示世界坐标系下的三维点.
		/* matR表示旋转矩阵.
		/* vecT表示平移向量.
		/* matQp 表示输出结果.
		/*******************************/
		int xFormProj(const Eigen::Matrix3Xd &matP, const Eigen::Matrix3d &matR,
			const Eigen::Vector3d &vecT, Eigen::Matrix2Xd &matQp);


		/********************************/
		/* t表示输出向量
		/* R 表示旋转矩阵.
		/* G是3x3矩阵.
		/* F是3x3 集合.
		/* P 是三维点.
		/*******************************/

		int estimate_t(Eigen::Vector3d &t, const Eigen::Matrix3d &R, const Eigen::Matrix3d &G,
			const std::vector<Eigen::Matrix3d> &F, const Eigen::Matrix3Xd &P);


		int absKernel(Eigen::Matrix3d &R, Eigen::Vector3d &t, Eigen::Matrix3Xd &Qout,
			double &err2, const Eigen::Matrix3Xd &P, const Eigen::Matrix3Xd &Q,
			const std::vector<Eigen::Matrix3d> &F, const Eigen::Matrix3d &G);



	};

}//namespace cvg;


#endif
