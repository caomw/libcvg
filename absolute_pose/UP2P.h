#pragma once

#include <Eigen/Eigen>

#include "types_up2p.h"

namespace cvg
{

	class CUP2P
	{


	public:

		CUP2P();

		~CUP2P();

	public:

		void ComputeRT(
			Eigen::Matrix3Xd matXYZ,
			Eigen::Matrix2Xd matXY,
			Eigen::Matrix3d& R,
			Eigen::Vector3d& t);

	private:
		int up2p(
			const std::vector<Eigen::Vector3d>& x,
			const std::vector<Eigen::Vector3d>& X,
			CameraPoseVector* output);

	};

}
