

#ifndef __CP4PF_HPP__
#define __CP4PF_HPP__

#pragma warning(disable:4996)

#include <vector>

#include <Eigen/Eigen>

#include "types.h"

namespace cvg
{

	class CP4Pf
	{
	public:
		CP4Pf();
		virtual ~CP4Pf();

	public:
		void ComputeRT(
			Eigen::Matrix3Xd matXYZ,
			Eigen::Matrix2Xd matXY,
			Eigen::Matrix3d &R,
			Eigen::Vector3d &t);


	private:

		// Solves for camera pose and focal length(alpha) such that : lambda * diag(1, 1, alpha)*x = R * X + t
			// Re-implementation of the p4pf solver from
			//    Kukelova et al., Efficient Intersection of Three Quadrics and Applications in Computer Vision, CVPR 2016
			// Note that this solver does not enforce that the rows of the rotation are consistent. This also be interpreted as
			// having non-unit aspect ratio, i.e. fx = f * R.row(0).norm() and fy = f * R.row(1).norm();
		int p4pf(
			const std::vector<Eigen::Vector3d> &X,
			const std::vector<Eigen::Vector3d> &x,
			std::vector<CameraPose> *output);



	};

}//namespace cvg;






#endif
