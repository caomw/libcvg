#pragma once

#include <Eigen/Eigen>


namespace cvg
{

	class CUmeyama
	{

	public:
		CUmeyama();

		~CUmeyama();


	public:

		/*********************************/
		/* matU 表示三维点.
		/* matv 表示二维点.
		/* poseR表示旋转矩阵.
		/* poseT 表示平移向量.
		/*********************************/
		void ComputeRT(
			Eigen::Matrix3Xd matU, Eigen::Matrix2Xd matv,
			Eigen::Matrix3d &poseR, Eigen::Vector3d &poseT);

	};

}// namespace cvg;
