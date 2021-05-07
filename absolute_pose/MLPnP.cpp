

#pragma warning(disable:4996)


#include "MLPnP.hpp"


#include <iostream>
#include <vector>



namespace cvg
{


	CMLPnP::CMLPnP()
	{


	}


	CMLPnP::~CMLPnP()
	{

	}

	int CMLPnP::ComputeRT(
		Eigen::Matrix3d& mR, Eigen::Vector3d& mT,
		Eigen::Matrix3Xd& MatXYZ, Eigen::Matrix2Xd& Matxy, double tol /*= 1e-8*/)
	{
		if (MatXYZ.cols() == 0 || MatXYZ.rows() == 0)
		{
			std::cout << "MatXYZ is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		if (Matxy.cols() == 0 || Matxy.rows() == 0)
		{
			std::cout << "Matxy is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}


		if (MatXYZ.cols() != Matxy.cols())
		{
			std::cout << "The columns of MatXYZ is not equal to that of Matxy !"
				<< " at line: " << __LINE__ << ", in file: " << __FILE__ << "." << std::endl;
			return -1;
		}


		std::size_t nrPts = Matxy.cols();

		Eigen::Matrix3Xd tempV = Eigen::Matrix3Xd::Ones(3, nrPts);
		tempV.topRows<2>() = Matxy;

		Eigen::RowVectorXd normV = tempV.colwise().norm();


		tempV.row(0) = tempV.row(0).cwiseQuotient(normV);
		tempV.row(1) = tempV.row(1).cwiseQuotient(normV);
		tempV.row(2) = tempV.row(2).cwiseQuotient(normV);

		Eigen::MatrixXd matT;
		Statistics statistics;


		int flag = MLPnPSolver_Without_Cov(MatXYZ, tempV, matT, statistics);


		if (-1 == flag)
		{
			std::cout << " Error happend when invoke MLPnPSolver function !"
				<< " at line: " << __LINE__ << " in file: " << __FILE__ << "." << std::endl;

			return -1;
		}


		mR = matT.leftCols<3>();
		mT = matT.rightCols<1>();


		return 1;

	}

	int CMLPnP::RodrigueToRotM(Eigen::Vector3d& omega, Eigen::Matrix3d& RotM)
	{
		if (omega.size() != 3)
		{
			std::cout << "The omega is empty, " << "at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;

		}

		RotM.setZero();

		Eigen::Matrix3d wx = Eigen::Matrix3d::Zero();

		wx(0, 1) = -omega(2);
		wx(0, 2) = omega(1);

		wx(1, 0) = omega(2);
		wx(1, 2) = -omega(0);

		wx(2, 0) = -omega(1);
		wx(2, 1) = omega(0);

		double omega_norm = omega.norm();

		double eps = std::numeric_limits<double>::epsilon();

		if (omega_norm < eps)
		{
			RotM = Eigen::Matrix3d::Identity();
		}
		else
		{

			RotM = Eigen::Matrix3d::Identity();
			RotM = RotM + std::sin(omega_norm) / omega_norm * wx + (1 - std::cos(omega_norm)) / (omega_norm * omega_norm) * (wx * wx);

		}


		return 1;
	}





	int CMLPnP::RodrigueToVect(Eigen::Matrix3d& RotM, Eigen::Vector3d& omega)
	{
		if (RotM.size() != 9)
		{
			std::cout << "The RotM is empty or incorrect " << ", at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;
		}


		double w_norm = std::acos((RotM.trace() - 1) / 2.0);

		double eps = std::numeric_limits<double>::epsilon();


		if (w_norm < eps)
		{
			omega = Eigen::Vector3d::Zero();
		}
		else
		{
			Eigen::Vector3d temp;
			temp(0) = RotM(2, 1) - RotM(1, 2);
			temp(1) = RotM(0, 2) - RotM(2, 0);
			temp(2) = RotM(1, 0) - RotM(0, 1);

			omega = 1.0 / (2 * std::sin(w_norm)) * temp * w_norm;

		}

		return 1;

	}

	int CMLPnP::JacobiansRodrigues(double X1, double Y1, double Z1,
		double r1, double r2, double r3, double s1, double s2, double s3,
		double t1, double t2, double t3, double w1, double w2, double w3, Eigen::MatrixXd& jacM)
	{


		double t5 = w1 * w1;
		double t6 = w2 * w2;
		double t7 = w3 * w3;
		double t8 = t5 + t6 + t7;
		double t9 = std::sqrt(t8);
		double t10 = std::sin(t9);
		double t11 = 1.0 / std::sqrt(t8);
		double t12 = std::cos(t9);
		double t13 = t12 - 1.0;
		double t14 = 1.0 / t8;
		double t16 = t10 * t11 * w3;
		double t17 = t13 * t14 * w1 * w2;
		double t19 = t10 * t11 * w2;
		double t20 = t13 * t14 * w1 * w3;
		double t24 = t6 + t17;
		double t27 = t16 + t17;
		double t28 = Y1 * t27;
		double t29 = t19 - t20;
		double t30 = Z1 * t29;
		double t31 = t13 * t14 * t24;
		double t32 = t31 + 1.0;
		double t33 = X1 * t32;
		double t15 = t1 - t28 + t30 + t33;
		double t21 = t10 * t11 * w1;
		double t22 = t13 * t14 * w2 * w3;
		double t45 = t5 + t7;
		double t53 = t16 - t17;
		double t54 = X1 * t53;
		double t55 = t21 + t22;
		double t56 = Z1 * t55;
		double t57 = t13 * t14 * t45;
		double t58 = t57 + 1.0;
		double t59 = Y1 * t58;
		double t18 = t2 + t54 - t56 + t59;
		double t34 = t5 + t6;
		double t38 = t19 + t20;
		double t39 = X1 * t38;
		double t40 = t21 - t22;
		double t41 = Y1 * t40;
		double t42 = t13 * t14 * t34;
		double t43 = t42 + 1.0;
		double t44 = Z1 * t43;
		double t23 = t3 - t39 + t41 + t44;
		double t25 = 1.0 / std::pow(t8, 3.0 / 2.0);
		double t26 = 1.0 / (t8 * t8);
		double t35 = t12 * t14 * w1 * w2;
		double t36 = t5 * t10 * t25 * w3;
		double t37 = t5 * t13 * t26 * w3 * 2.0;
		double t46 = t10 * t25 * w1 * w3;
		double t47 = t5 * t10 * t25 * w2;
		double t48 = t5 * t13 * t26 * w2 * 2.0;
		double t49 = t10 * t11;
		double t50 = t5 * t12 * t14;
		double t51 = t13 * t26 * w1 * w2 * w3 * 2.0;
		double t52 = t10 * t25 * w1 * w2 * w3;
		double t60 = t15 * t15;
		double t61 = t18 * t18;
		double t62 = t23 * t23;
		double t63 = t60 + t61 + t62;
		double t64 = t5 * t10 * t25;
		double t65 = 1.0 / std::sqrt(t63);
		double t66 = Y1 * r2 * t6;
		double t67 = Z1 * r3 * t7;
		double t68 = r1 * t1 * t5;
		double t69 = r1 * t1 * t6;
		double t70 = r1 * t1 * t7;
		double t71 = r2 * t2 * t5;
		double t72 = r2 * t2 * t6;
		double t73 = r2 * t2 * t7;
		double t74 = r3 * t3 * t5;
		double t75 = r3 * t3 * t6;
		double t76 = r3 * t3 * t7;
		double t77 = X1 * r1 * t5;
		double t78 = X1 * r2 * w1 * w2;
		double t79 = X1 * r3 * w1 * w3;
		double t80 = Y1 * r1 * w1 * w2;
		double t81 = Y1 * r3 * w2 * w3;
		double t82 = Z1 * r1 * w1 * w3;
		double t83 = Z1 * r2 * w2 * w3;
		double t84 = X1 * r1 * t6 * t12;
		double t85 = X1 * r1 * t7 * t12;
		double t86 = Y1 * r2 * t5 * t12;
		double t87 = Y1 * r2 * t7 * t12;
		double t88 = Z1 * r3 * t5 * t12;
		double t89 = Z1 * r3 * t6 * t12;
		double t90 = X1 * r2 * t9 * t10 * w3;
		double t91 = Y1 * r3 * t9 * t10 * w1;
		double t92 = Z1 * r1 * t9 * t10 * w2;
		double t102 = X1 * r3 * t9 * t10 * w2;
		double t103 = Y1 * r1 * t9 * t10 * w3;
		double t104 = Z1 * r2 * t9 * t10 * w1;
		double t105 = X1 * r2 * t12 * w1 * w2;
		double t106 = X1 * r3 * t12 * w1 * w3;
		double t107 = Y1 * r1 * t12 * w1 * w2;
		double t108 = Y1 * r3 * t12 * w2 * w3;
		double t109 = Z1 * r1 * t12 * w1 * w3;
		double t110 = Z1 * r2 * t12 * w2 * w3;

		double t93 = t66 + t67 + t68 + t69 + t70 + t71 + t72 + t73 + t74 + t75 + t76 + t77 + t78 + t79 + t80 + t81 + t82 +
			t83 + t84 + t85 + t86 + t87 + t88 + t89 + t90 + t91 + t92 - t102 - t103 - t104 - t105 - t106 - t107 - t108 - t109 - t110;

		double t94 = t10 * t25 * w1 * w2;
		double t95 = t6 * t10 * t25 * w3;
		double t96 = t6 * t13 * t26 * w3 * 2.0;
		double t97 = t12 * t14 * w2 * w3;
		double t98 = t6 * t10 * t25 * w1;
		double t99 = t6 * t13 * t26 * w1 * 2.0;
		double t100 = t6 * t10 * t25;
		double t101 = 1.0 / std::pow(t63, (3.0 / 2.0));
		double t111 = t6 * t12 * t14;
		double t112 = t10 * t25 * w2 * w3;
		double t113 = t12 * t14 * w1 * w3;
		double t114 = t7 * t10 * t25 * w2;
		double t115 = t7 * t13 * t26 * w2 * 2.0;
		double t116 = t7 * t10 * t25 * w1;
		double t117 = t7 * t13 * t26 * w1 * 2.0;
		double t118 = t7 * t12 * t14;
		double t119 = t13 * t24 * t26 * w1 * 2.0;
		double t120 = t10 * t24 * t25 * w1;

		double t121 = t119 + t120;
		double t122 = t13 * t26 * t34 * w1 * 2.0;
		double t123 = t10 * t25 * t34 * w1;
		double t131 = t13 * t14 * w1 * 2.0;
		double t124 = t122 + t123 - t131;
		double t139 = t13 * t14 * w3;
		double t125 = -t35 + t36 + t37 + t94 - t139;
		double t126 = X1 * t125;
		double t127 = t49 + t50 + t51 + t52 - t64;
		double t128 = Y1 * t127;
		double t129 = t126 + t128 - Z1 * t124;
		double t130 = t23 * t129 * 2.0;
		double t132 = t13 * t26 * t45 * w1 * 2.0;
		double t133 = t10 * t25 * t45 * w1;
		double t138 = t13 * t14 * w2;
		double t134 = -t46 + t47 + t48 + t113 - t138;
		double t135 = X1 * t134;
		double t136 = -t49 - t50 + t51 + t52 + t64;
		double t137 = Z1 * t136;
		double t140 = X1 * s1 * t5;
		double t141 = Y1 * s2 * t6;
		double t142 = Z1 * s3 * t7;
		double t143 = s1 * t1 * t5;
		double t144 = s1 * t1 * t6;
		double t145 = s1 * t1 * t7;
		double t146 = s2 * t2 * t5;
		double t147 = s2 * t2 * t6;
		double t148 = s2 * t2 * t7;
		double t149 = s3 * t3 * t5;
		double t150 = s3 * t3 * t6;
		double t151 = s3 * t3 * t7;
		double t152 = X1 * s2 * w1 * w2;
		double t153 = X1 * s3 * w1 * w3;
		double t154 = Y1 * s1 * w1 * w2;
		double t155 = Y1 * s3 * w2 * w3;
		double t156 = Z1 * s1 * w1 * w3;
		double t157 = Z1 * s2 * w2 * w3;
		double t158 = X1 * s1 * t6 * t12;
		double t159 = X1 * s1 * t7 * t12;
		double t160 = Y1 * s2 * t5 * t12;
		double t161 = Y1 * s2 * t7 * t12;
		double t162 = Z1 * s3 * t5 * t12;
		double t163 = Z1 * s3 * t6 * t12;
		double t164 = X1 * s2 * t9 * t10 * w3;
		double t165 = Y1 * s3 * t9 * t10 * w1;
		double t166 = Z1 * s1 * t9 * t10 * w2;
		double t183 = X1 * s3 * t9 * t10 * w2;
		double t184 = Y1 * s1 * t9 * t10 * w3;
		double t185 = Z1 * s2 * t9 * t10 * w1;
		double t186 = X1 * s2 * t12 * w1 * w2;
		double t187 = X1 * s3 * t12 * w1 * w3;
		double t188 = Y1 * s1 * t12 * w1 * w2;
		double t189 = Y1 * s3 * t12 * w2 * w3;
		double t190 = Z1 * s1 * t12 * w1 * w3;
		double t191 = Z1 * s2 * t12 * w2 * w3;

		double t167 = t140 + t141 + t142 + t143 + t144 + t145 + t146 + t147 + t148 + t149 + t150 + t151 + t152 +
			t153 + t154 + t155 + t156 + t157 + t158 + t159 + t160 + t161 + t162 + t163 + t164 + t165 + t166 -
			t183 - t184 - t185 - t186 - t187 - t188 - t189 - t190 - t191;

		double t168 = t13 * t26 * t45 * w2 * 2.0;
		double t169 = t10 * t25 * t45 * w2;
		double t170 = t168 + t169;
		double t171 = t13 * t26 * t34 * w2 * 2.0;
		double t172 = t10 * t25 * t34 * w2;
		double t176 = t13 * t14 * w2 * 2.0;
		double t173 = t171 + t172 - t176;
		double t174 = -t49 + t51 + t52 + t100 - t111;
		double t175 = X1 * t174;
		double t177 = t13 * t24 * t26 * w2 * 2.0;
		double t178 = t10 * t24 * t25 * w2;
		double t192 = t13 * t14 * w1;
		double t179 = -t97 + t98 + t99 + t112 - t192;
		double t180 = Y1 * t179;
		double t181 = t49 + t51 + t52 - t100 + t111;
		double t182 = Z1 * t181;
		double t193 = t13 * t26 * t34 * w3 * 2.0;
		double t194 = t10 * t25 * t34 * w3;
		double t195 = t193 + t194;
		double t196 = t13 * t26 * t45 * w3 * 2.0;
		double t197 = t10 * t25 * t45 * w3;
		double t200 = t13 * t14 * w3 * 2.0;
		double t198 = t196 + t197 - t200;
		double t199 = t7 * t10 * t25;
		double t201 = t13 * t24 * t26 * w3 * 2.0;
		double t202 = t10 * t24 * t25 * w3;
		double t203 = -t49 + t51 + t52 - t118 + t199;
		double t204 = Y1 * t203;
		double t205 = t1 * 2.0;
		double t206 = Z1 * t29 * 2.0;
		double t207 = X1 * t32 * 2.0;
		double t208 = t205 + t206 + t207 - Y1 * t27 * 2.0;
		double t209 = t2 * 2.0;
		double t210 = X1 * t53 * 2.0;
		double t211 = Y1 * t58 * 2.0;
		double t212 = t209 + t210 + t211 - Z1 * t55 * 2.0;
		double t213 = t3 * 2.0;
		double t214 = Y1 * t40 * 2.0;
		double t215 = Z1 * t43 * 2.0;
		double t216 = t213 + t214 + t215 - X1 * t38 * 2.0;


		double temp1 = t14 * t65 * (X1 * r1 * w1 * 2.0 + X1 * r2 * w2 + X1 * r3 * w3 + Y1 * r1 * w2 +
			Z1 * r1 * w3 + r1 * t1 * w1 * 2.0 + r2 * t2 * w1 * 2.0 + r3 * t3 * w1 * 2.0 + Y1 * r3 * t5 * t12 +
			Y1 * r3 * t9 * t10 - Z1 * r2 * t5 * t12 - Z1 * r2 * t9 * t10 - X1 * r2 * t12 * w2 - X1 * r3 * t12 * w3 -
			Y1 * r1 * t12 * w2 + Y1 * r2 * t12 * w1 * 2.0 - Z1 * r1 * t12 * w3 + Z1 * r3 * t12 * w1 * 2.0 +
			Y1 * r3 * t5 * t10 * t11 - Z1 * r2 * t5 * t10 * t11 + X1 * r2 * t12 * w1 * w3 -
			X1 * r3 * t12 * w1 * w2 - Y1 * r1 * t12 * w1 * w3 + Z1 * r1 * t12 * w1 * w2 -
			Y1 * r1 * t10 * t11 * w1 * w3 + Z1 * r1 * t10 * t11 * w1 * w2 -
			X1 * r1 * t6 * t10 * t11 * w1 - X1 * r1 * t7 * t10 * t11 * w1 + X1 * r2 * t5 * t10 * t11 * w2 +
			X1 * r3 * t5 * t10 * t11 * w3 + Y1 * r1 * t5 * t10 * t11 * w2 - Y1 * r2 * t5 * t10 * t11 * w1 -
			Y1 * r2 * t7 * t10 * t11 * w1 + Z1 * r1 * t5 * t10 * t11 * w3 - Z1 * r3 * t5 * t10 * t11 * w1 -
			Z1 * r3 * t6 * t10 * t11 * w1 + X1 * r2 * t10 * t11 * w1 * w3 - X1 * r3 * t10 * t11 * w1 * w2 +
			Y1 * r3 * t10 * t11 * w1 * w2 * w3 + Z1 * r2 * t10 * t11 * w1 * w2 * w3) - t26 * t65 * t93 * w1 * 2.0 -
			t14 * t93 * t101 * (t130 + t15 * (-X1 * t121 + Y1 * (t46 + t47 + t48 - t13 * t14 * w2 - t12 * t14 * w1 * w3) +
				Z1 * (t35 + t36 + t37 - t13 * t14 * w3 - t10 * t25 * w1 * w2)) * 2.0 +
				t18 * (t135 + t137 - Y1 * (t132 + t133 - t13 * t14 * w1 * 2.0)) * 2.0) * (1.0 / 2.0);

		double temp2 = t14 * t65 * (X1 * s1 * w1 * 2.0 + X1 * s2 * w2 +
			X1 * s3 * w3 + Y1 * s1 * w2 + Z1 * s1 * w3 + s1 * t1 * w1 * 2.0 + s2 * t2 * w1 * 2.0 +
			s3 * t3 * w1 * 2.0 + Y1 * s3 * t5 * t12 + Y1 * s3 * t9 * t10 - Z1 * s2 * t5 * t12 - Z1 * s2 * t9 * t10 - X1 * s2 * t12 * w2 -
			X1 * s3 * t12 * w3 - Y1 * s1 * t12 * w2 + Y1 * s2 * t12 * w1 * 2.0 - Z1 * s1 * t12 * w3 +
			Z1 * s3 * t12 * w1 * 2.0 + Y1 * s3 * t5 * t10 * t11 - Z1 * s2 * t5 * t10 * t11 +
			X1 * s2 * t12 * w1 * w3 - X1 * s3 * t12 * w1 * w2 - Y1 * s1 * t12 * w1 * w3 +
			Z1 * s1 * t12 * w1 * w2 + X1 * s2 * t10 * t11 * w1 * w3 - X1 * s3 * t10 * t11 * w1 * w2 -
			Y1 * s1 * t10 * t11 * w1 * w3 + Z1 * s1 * t10 * t11 * w1 * w2 - X1 * s1 * t6 * t10 * t11 * w1 -
			X1 * s1 * t7 * t10 * t11 * w1 + X1 * s2 * t5 * t10 * t11 * w2 + X1 * s3 * t5 * t10 * t11 * w3 +
			Y1 * s1 * t5 * t10 * t11 * w2 - Y1 * s2 * t5 * t10 * t11 * w1 - Y1 * s2 * t7 * t10 * t11 * w1 +
			Z1 * s1 * t5 * t10 * t11 * w3 - Z1 * s3 * t5 * t10 * t11 * w1 - Z1 * s3 * t6 * t10 * t11 * w1 +
			Y1 * s3 * t10 * t11 * w1 * w2 * w3 + Z1 * s2 * t10 * t11 * w1 * w2 * w3) -
			t14 * t101 * t167 * (t130 + t15 * (Y1 * (t46 + t47 + t48 - t113 - t138) +
				Z1 * (t35 + t36 + t37 - t94 - t139) - X1 * t121) * 2.0 + t18 * (t135 + t137 -
					Y1 * (-t131 + t132 + t133)) * 2.0) * (1.0 / 2.0) - t26 * t65 * t167 * w1 * 2.0;

		double temp3 = t14 * t65 * (X1 * r2 * w1 +
			Y1 * r1 * w1 + Y1 * r2 * w2 * 2.0 + Y1 * r3 * w3 + Z1 * r2 * w3 + r1 * t1 * w2 * 2.0 +
			r2 * t2 * w2 * 2.0 + r3 * t3 * w2 * 2.0 - X1 * r3 * t6 * t12 - X1 * r3 * t9 * t10 +
			Z1 * r1 * t6 * t12 + Z1 * r1 * t9 * t10 + X1 * r1 * t12 * w2 * 2.0 - X1 * r2 * t12 * w1 -
			Y1 * r1 * t12 * w1 - Y1 * r3 * t12 * w3 - Z1 * r2 * t12 * w3 + Z1 * r3 * t12 * w2 * 2.0 -
			X1 * r3 * t6 * t10 * t11 + Z1 * r1 * t6 * t10 * t11 + X1 * r2 * t12 * w2 * w3 - Y1 * r1 * t12 * w2 * w3 +
			Y1 * r3 * t12 * w1 * w2 - Z1 * r2 * t12 * w1 * w2 - Y1 * r1 * t10 * t11 * w2 * w3 +
			Y1 * r3 * t10 * t11 * w1 * w2 - Z1 * r2 * t10 * t11 * w1 * w2 - X1 * r1 * t6 * t10 * t11 * w2 +
			X1 * r2 * t6 * t10 * t11 * w1 - X1 * r1 * t7 * t10 * t11 * w2 + Y1 * r1 * t6 * t10 * t11 * w1 -
			Y1 * r2 * t5 * t10 * t11 * w2 - Y1 * r2 * t7 * t10 * t11 * w2 + Y1 * r3 * t6 * t10 * t11 * w3 -
			Z1 * r3 * t5 * t10 * t11 * w2 + Z1 * r2 * t6 * t10 * t11 * w3 - Z1 * r3 * t6 * t10 * t11 * w2 +
			X1 * r2 * t10 * t11 * w2 * w3 + X1 * r3 * t10 * t11 * w1 * w2 * w3 + Z1 * r1 * t10 * t11 * w1 * w2 * w3) -
			t26 * t65 * t93 * w2 * 2.0 - t14 * t93 * t101 * (t18 * (Z1 * (-t35 + t94 + t95 + t96 - t13 * t14 * w3) -
				Y1 * t170 + X1 * (t97 + t98 + t99 - t13 * t14 * w1 - t10 * t25 * w2 * w3)) * 2.0 +
				t15 * (t180 + t182 - X1 * (t177 + t178 - t13 * t14 * w2 * 2.0)) * 2.0 + t23 * (t175 +
					Y1 * (t35 - t94 + t95 + t96 - t13 * t14 * w3) - Z1 * t173) * 2.0) * (1.0 / 2.0);

		double temp4 = t14 * t65 * (X1 * s2 * w1 +
			Y1 * s1 * w1 + Y1 * s2 * w2 * 2.0 + Y1 * s3 * w3 + Z1 * s2 * w3 + s1 * t1 * w2 * 2.0 + s2 * t2 * w2 * 2.0 +
			s3 * t3 * w2 * 2.0 - X1 * s3 * t6 * t12 - X1 * s3 * t9 * t10 + Z1 * s1 * t6 * t12 + Z1 * s1 * t9 * t10 +
			X1 * s1 * t12 * w2 * 2.0 - X1 * s2 * t12 * w1 - Y1 * s1 * t12 * w1 - Y1 * s3 * t12 * w3 - Z1 * s2 * t12 * w3 +
			Z1 * s3 * t12 * w2 * 2.0 - X1 * s3 * t6 * t10 * t11 + Z1 * s1 * t6 * t10 * t11 + X1 * s2 * t12 * w2 * w3 -
			Y1 * s1 * t12 * w2 * w3 + Y1 * s3 * t12 * w1 * w2 - Z1 * s2 * t12 * w1 * w2 + X1 * s2 * t10 * t11 * w2 * w3 -
			Y1 * s1 * t10 * t11 * w2 * w3 + Y1 * s3 * t10 * t11 * w1 * w2 - Z1 * s2 * t10 * t11 * w1 * w2 -
			X1 * s1 * t6 * t10 * t11 * w2 + X1 * s2 * t6 * t10 * t11 * w1 - X1 * s1 * t7 * t10 * t11 * w2 +
			Y1 * s1 * t6 * t10 * t11 * w1 - Y1 * s2 * t5 * t10 * t11 * w2 - Y1 * s2 * t7 * t10 * t11 * w2 +
			Y1 * s3 * t6 * t10 * t11 * w3 - Z1 * s3 * t5 * t10 * t11 * w2 + Z1 * s2 * t6 * t10 * t11 * w3 -
			Z1 * s3 * t6 * t10 * t11 * w2 + X1 * s3 * t10 * t11 * w1 * w2 * w3 + Z1 * s1 * t10 * t11 * w1 * w2 * w3) -
			t26 * t65 * t167 * w2 * 2.0 - t14 * t101 * t167 * (t18 * (X1 * (t97 + t98 + t99 - t112 - t192) +
				Z1 * (-t35 + t94 + t95 + t96 - t139) - Y1 * t170) * 2.0 + t15 * (t180 + t182 - X1 * (-t176 + t177 + t178)) * 2.0 +
				t23 * (t175 + Y1 * (t35 - t94 + t95 + t96 - t139) - Z1 * t173) * 2.0) * (1.0 / 2.0);


		double temp5 = t14 * t65 * (X1 * r3 * w1 +
			Y1 * r3 * w2 + Z1 * r1 * w1 + Z1 * r2 * w2 + Z1 * r3 * w3 * 2.0 + r1 * t1 * w3 * 2.0 + r2 * t2 * w3 * 2.0 +
			r3 * t3 * w3 * 2.0 + X1 * r2 * t7 * t12 + X1 * r2 * t9 * t10 - Y1 * r1 * t7 * t12 - Y1 * r1 * t9 * t10 +
			X1 * r1 * t12 * w3 * 2.0 - X1 * r3 * t12 * w1 + Y1 * r2 * t12 * w3 * 2.0 - Y1 * r3 * t12 * w2 -
			Z1 * r1 * t12 * w1 - Z1 * r2 * t12 * w2 + X1 * r2 * t7 * t10 * t11 - Y1 * r1 * t7 * t10 * t11 -
			X1 * r3 * t12 * w2 * w3 + Y1 * r3 * t12 * w1 * w3 + Z1 * r1 * t12 * w2 * w3 - Z1 * r2 * t12 * w1 * w3 +
			Y1 * r3 * t10 * t11 * w1 * w3 + Z1 * r1 * t10 * t11 * w2 * w3 - Z1 * r2 * t10 * t11 * w1 * w3 -
			X1 * r1 * t6 * t10 * t11 * w3 - X1 * r1 * t7 * t10 * t11 * w3 + X1 * r3 * t7 * t10 * t11 * w1 -
			Y1 * r2 * t5 * t10 * t11 * w3 - Y1 * r2 * t7 * t10 * t11 * w3 + Y1 * r3 * t7 * t10 * t11 * w2 +
			Z1 * r1 * t7 * t10 * t11 * w1 + Z1 * r2 * t7 * t10 * t11 * w2 - Z1 * r3 * t5 * t10 * t11 * w3 -
			Z1 * r3 * t6 * t10 * t11 * w3 - X1 * r3 * t10 * t11 * w2 * w3 + X1 * r2 * t10 * t11 * w1 * w2 * w3 +
			Y1 * r1 * t10 * t11 * w1 * w2 * w3) - t26 * t65 * t93 * w3 * 2.0 - t14 * t93 * t101 * (t18 * (Z1 * (t46 -
				t113 + t114 + t115 - t13 * t14 * w2) - Y1 * t198 + X1 * (t49 + t51 + t52 + t118 - t7 * t10 * t25)) * 2.0 +
				t23 * (X1 * (-t97 + t112 + t116 + t117 - t13 * t14 * w1) + Y1 * (-t46 + t113 + t114 + t115 - t13 * t14 * w2) -
					Z1 * t195) * 2.0 + t15 * (t204 + Z1 * (t97 - t112 + t116 + t117 - t13 * t14 * w1) -
						X1 * (t201 + t202 - t13 * t14 * w3 * 2.0)) * 2.0) * (1.0 / 2.0);

		double temp6 = t14 * t65 * (X1 * s3 * w1 +
			Y1 * s3 * w2 + Z1 * s1 * w1 + Z1 * s2 * w2 + Z1 * s3 * w3 * 2.0 + s1 * t1 * w3 * 2.0 + s2 * t2 * w3 * 2.0 +
			s3 * t3 * w3 * 2.0 + X1 * s2 * t7 * t12 + X1 * s2 * t9 * t10 - Y1 * s1 * t7 * t12 - Y1 * s1 * t9 * t10 +
			X1 * s1 * t12 * w3 * 2.0 - X1 * s3 * t12 * w1 + Y1 * s2 * t12 * w3 * 2.0 -
			Y1 * s3 * t12 * w2 - Z1 * s1 * t12 * w1 - Z1 * s2 * t12 * w2 + X1 * s2 * t7 * t10 * t11 -
			Y1 * s1 * t7 * t10 * t11 - X1 * s3 * t12 * w2 * w3 + Y1 * s3 * t12 * w1 * w3 + Z1 * s1 * t12 * w2 * w3 -
			Z1 * s2 * t12 * w1 * w3 - X1 * s3 * t10 * t11 * w2 * w3 + Y1 * s3 * t10 * t11 * w1 * w3 +
			Z1 * s1 * t10 * t11 * w2 * w3 - Z1 * s2 * t10 * t11 * w1 * w3 - X1 * s1 * t6 * t10 * t11 * w3 -
			X1 * s1 * t7 * t10 * t11 * w3 + X1 * s3 * t7 * t10 * t11 * w1 - Y1 * s2 * t5 * t10 * t11 * w3 -
			Y1 * s2 * t7 * t10 * t11 * w3 + Y1 * s3 * t7 * t10 * t11 * w2 + Z1 * s1 * t7 * t10 * t11 * w1 +
			Z1 * s2 * t7 * t10 * t11 * w2 - Z1 * s3 * t5 * t10 * t11 * w3 - Z1 * s3 * t6 * t10 * t11 * w3 +
			X1 * s2 * t10 * t11 * w1 * w2 * w3 + Y1 * s1 * t10 * t11 * w1 * w2 * w3) - t26 * t65 * t167 * w3 * 2.0 -
			t14 * t101 * t167 * (t18 * (Z1 * (t46 - t113 + t114 + t115 - t138) - Y1 * t198 +
				X1 * (t49 + t51 + t52 + t118 - t199)) * 2.0 + t23 * (X1 * (-t97 + t112 + t116 + t117 -
					t192) + Y1 * (-t46 + t113 + t114 + t115 - t138) - Z1 * t195) * 2.0 + t15 * (t204 + Z1 * (t97 - t112 +
						t116 + t117 - t192) - X1 * (-t200 + t201 + t202)) * 2.0) * (1.0 / 2.0);

		double temp7 = r1 * t65 - t14 * t93 * t101 * t208 * (1.0 / 2.0);

		double temp8 = s1 * t65 - t14 * t101 * t167 * t208 * (1.0 / 2.0);

		double temp9 = r2 * t65 - t14 * t93 * t101 * t212 * (1.0 / 2.0);
		double temp10 = s2 * t65 - t14 * t101 * t167 * t212 * (1.0 / 2.0);
		double temp11 = r3 * t65 - t14 * t93 * t101 * t216 * (1.0 / 2.0);
		double temp12 = s3 * t65 - t14 * t101 * t167 * t216 * (1.0 / 2.0);


		jacM = Eigen::MatrixXd::Zero(2, 6);

		//第一列.
		jacM(0, 0) = temp1;
		jacM(1, 0) = temp2;

		//第二列.
		jacM(0, 1) = temp3;
		jacM(1, 1) = temp4;

		//第3列.
		jacM(0, 2) = temp5;
		jacM(1, 2) = temp6;

		//第4列.
		jacM(0, 3) = temp7;
		jacM(1, 3) = temp8;

		//第5列.
		jacM(0, 4) = temp9;
		jacM(1, 4) = temp10;

		//第6列.
		jacM(0, 5) = temp11;
		jacM(1, 5) = temp12;



		return 1;
	}

	int CMLPnP::ResidualsAndJacobian(Eigen::VectorXd& x, Eigen::Matrix3Xd& r,
		Eigen::Matrix3Xd& s, Eigen::Matrix3Xd& points3D, Eigen::VectorXd& err, Eigen::MatrixXd& J)
	{

		if (x.size() != 6)
		{
			std::cout << "x is incorrect ! " << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;
			return -1;
		}

		if (r.cols() <= 0)
		{
			std::cout << " r is incorrect !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;
			return -1;

		}

		if (s.cols() <= 0)
		{
			std::cout << " s is incorrect !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;
			return -1;
		}

		if (points3D.cols() <= 0)
		{
			std::cout << " points3D is incorrect !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;
			return -1;
		}

		std::size_t nrPts = points3D.cols();

		err = Eigen::VectorXd::Zero(2 * nrPts);

		J = Eigen::MatrixXd::Zero(2 * nrPts, 6);

		/////////////////////////////
		Eigen::Matrix3d R;
		Eigen::Vector3d omega;
		omega(0) = x(0);
		omega(1) = x(1);
		omega(2) = x(2);

		int flag1 = RodrigueToRotM(omega, R);

		if (-1 == flag1)
		{
			std::cout << "Error happend when invoke RodrigueToRotM function !"
				<< " at line: " << __LINE__ << ", in file: " << __FILE__ << "." << std::endl;
			return -1;
		}

		///////////////////////////////

		Eigen::Vector3d t = x.tail<3>();

		Eigen::Matrix3Xd temp = Eigen::Matrix3Xd::Zero(3, nrPts);
		temp.colwise() = t;

		Eigen::Matrix3Xd res1 = R * points3D + temp;

		Eigen::RowVectorXd normc = res1.colwise().norm();

		Eigen::Matrix3Xd normres = Eigen::Matrix3Xd::Zero(3, nrPts);
		normres.row(0) = res1.row(0).cwiseQuotient(normc);
		normres.row(1) = res1.row(1).cwiseQuotient(normc);
		normres.row(2) = res1.row(2).cwiseQuotient(normc);


		err = Eigen::VectorXd::Zero(2 * nrPts);

		for (std::size_t i = 0; i < r.cols();/*nrPts*/ i++)
		{
			err(2 * i) = r.col(i).transpose() * normres.col(i);
			err(2 * i + 1) = s.col(i).transpose() * normres.col(i);

			Eigen::MatrixXd jacobian;
			int flag = JacobiansRodrigues(points3D(0, i), points3D(1, i), points3D(2, i),
				r(0, i), r(1, i), r(2, i), s(0, i), s(1, i), s(2, i), x(3), x(4), x(5), x(0), x(1), x(2), jacobian);

			if (flag == -1)
			{
				std::cout << "error happend !" << " at line: "
					<< __LINE__ << ", in file: " << __FILE__ << std::endl;
				return -1;
			}

			J.row(2 * i) = jacobian.row(0);
			J.row(2 * i + 1) = jacobian.row(1);

		}

		return 1;
	}

	int CMLPnP::OptimMLPnP_GN(
		Eigen::Matrix<double, 3, 4, Eigen::ColMajor>& Tinit,
		Eigen::Matrix3Xd& points3D, Eigen::Matrix3Xd& rnull,
		Eigen::Matrix3Xd& snull, Eigen::MatrixXd& P,
		const OptimFlags& optimFlags,
		Eigen::Matrix<double, 3, 4, Eigen::ColMajor>& Tout, Statistics& statistics)
	{
		std::size_t cnt1 = points3D.cols();
		std::size_t cnt2 = rnull.cols();
		std::size_t cnt3 = snull.cols();
		std::size_t cnt4 = P.cols();

		if (cnt1 != cnt2 || cnt2 != cnt3 || cnt3 * 2 != cnt4)
		{
			std::cout << "The input parameters are incorrect !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;

		}

		// homogeneous to minimal
		Eigen::Matrix3d temp1 = Tinit.leftCols<3>();

		Eigen::VectorXd x = Eigen::VectorXd::Zero(6);


		Eigen::Vector3d omega;

		int flag = RodrigueToVect(temp1, omega);

		if (flag == -1)
		{
			std::cout << " Error happend when invoke RodrigueToVect function !"
				<< " at line: " << __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}

		x.head<3>() = omega;
		x.tail<3>() = Tinit.col(3);

		std::size_t nrl = rnull.cols();

		std::size_t redundanz = 2 * nrl - x.size();

		// optim params
		double epsParam = optimFlags.epsP;
		double epsFunc = optimFlags.epsF;

		// iteration params

		std::size_t cnt = 0;
		bool stop = false;
		Eigen::MatrixXd invKll = P;

		Eigen::VectorXd r;
		Eigen::MatrixXd J;

		Eigen::Matrix<double, 6, 6, Eigen::ColMajor> N;

		while (cnt < optimFlags.maxit && stop == 0)
		{

			int flag = ResidualsAndJacobian(x, rnull, snull, points3D, r, J);

			// design matrix
			N = J.transpose() * invKll * J;

			// System matrix
			Eigen::VectorXd g = J.transpose() * invKll * r;

			Eigen::VectorXd dx = N.inverse() * g;

			if (dx.cwiseAbs().maxCoeff() > 20 || dx.cwiseAbs().minCoeff() > 1)
			{
				break;
			}

			Eigen::VectorXd dl = J * dx;

			if (dl.cwiseAbs().maxCoeff() < epsFunc || dx.cwiseAbs().maxCoeff() < epsParam)
			{
				x = x - dx;
				break;
			}
			else
			{
				//update parameter vector

				x = x - dx;
			}

			cnt++;

		}// end while.

		// minimal to homogeneous

		Eigen::Vector3d omegaV = x.head<3>();
		Eigen::Matrix3d rotM;
		int flag2 = RodrigueToRotM(omegaV, rotM);

		if (flag2 == -1)
		{
			std::cout << "Error happend when invoke RodrigueToRotM function !"
				<< " at line: " << __LINE__ << " in file: " << __FILE__ << std::endl;
			return -1;
		}


		Tout.leftCols<3>() = rotM;
		Tout.col(3) = x.tail<3>();

		// empirical variance factor

		double resV = r.transpose() * invKll * r;

		double s0 = 0.0;

		if (redundanz > 0)
		{
			if (redundanz < nrl)
			{
				s0 = 1;
			}
			else
			{
				s0 = resV / redundanz;
			}


		}
		else
		{
			s0 = std::numeric_limits<double>::epsilon();
		}

		// variance - covariance matrix

		Eigen::Matrix<double, 6, 6, Eigen::ColMajor> Qxx = N.inverse();

		// cofactor matrix of "adjusted observations"
		Eigen::MatrixXd Qldld = J * statistics.Qxx * J.transpose();

		statistics.resV = resV;
		statistics.r = r;
		statistics.J = J;
		statistics.Qxx = Qxx;
		statistics.s0 = s0;
		statistics.Qldld = Qldld;
		statistics.param = (s0 * Qxx.diagonal()).cwiseSqrt();

		return 1;

	}// end function (OptimMLPnP_GN)



	/******************************************************/
	/* MLPnP solver.
	/* points3D, 输入参数, 表示世界坐标系下的三维点.
	/* v, 输入参数, 图像点(齐次坐标表示)
	/* T, 输出参数, 表示3x4矩阵.
	/* statistics, 输出参数,
	/*******************************************************/

	int CMLPnP::MLPnPSolver_Without_Cov(
		Eigen::Matrix3Xd& points3D, Eigen::Matrix3Xd& v,
		Eigen::MatrixXd& T, Statistics& statistics)
	{


		std::size_t nrPts = points3D.cols();
		std::size_t vn = v.cols();

		if (nrPts == 0)
		{
			std::cout << "points3D is empty !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;
		}

		if (vn == 0)
		{
			std::cout << "v is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}

		if (nrPts != vn)
		{
			std::cout << "The size of points3D is not same as that of v !"
				<< " at line: " << __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;

		}



		//matrix of null space vectors r and s

		Eigen::Matrix3Xd r = Eigen::Matrix3Xd::Zero(3, nrPts);
		Eigen::Matrix3Xd s = Eigen::Matrix3Xd::Zero(3, nrPts);

		std::vector<Eigen::Matrix2d> cov_reduced(nrPts);

		// test planarity, only works well if the scene is really planar
		// quasi - planar won't work very well

		Eigen::Matrix3d S = points3D * points3D.transpose();

		Eigen::EigenSolver<Eigen::Matrix3d> solver(S);

		Eigen::MatrixXd eigRot = solver.pseudoEigenvectors();


		int planar = 0;

		//create full design matrix

		Eigen::MatrixXd points3Dn;

		Eigen::MatrixXd A;

		//注意矩阵秩的求法.
		Eigen::FullPivLU <Eigen::Matrix3d> lu_decomp(S);

		if (lu_decomp.rank() == 2)
		{
			planar = 1;

			Eigen::Matrix2Xd points3D1 = eigRot.transpose() * points3D;

			points3Dn = Eigen::MatrixXd::Ones(3, nrPts);

			points3Dn.topRows<2>() = points3D1;

			//points3Dn 的第三行元素全为1.

			// create reduced design matrix

			A = Eigen::MatrixXd::Zero(2 * nrPts, 9);


		}
		else
		{

			points3Dn = Eigen::MatrixXd::Ones(4, nrPts);

			points3Dn.topRows<3>() = points3D;

			A = Eigen::MatrixXd::Zero(2 * nrPts, 12);

		}


		// compute null spaces of bearing vector v : null(v')

		//Eigen::MatrixXd null_2d = Eigen::MatrixXd::Zero(nrPts, 2);

		Eigen::MatrixXd null_2d;

		for (int i = 0; i < nrPts; i++)
		{

			//注意向量核空间的求解方法.
			// nullspace of right vector.
			Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner>
				svd_solver(v.col(i).transpose(), Eigen::ComputeFullV);

			//null_2d = svd_solver.matrixV().block(0, 1, 3, 2);
			null_2d = svd_solver.matrixV().rightCols<2>();

			r.col(i) = null_2d.col(0);
			s.col(i) = null_2d.col(1);

		}

		// stochastic model

		Eigen::MatrixXd Kll = Eigen::MatrixXd::Identity(2 * nrPts, 2 * nrPts);

		if (planar == 1)
		{
			// build reduces system.

			for (std::size_t i = 0; i < nrPts; i++)
			{
				//r12.			
				A(2 * i, 0) = r(0, i) * points3Dn(1, i);
				A(2 * i + 1, 0) = s(0, i) * points3Dn(1, i);

				//r13
				A(2 * i, 1) = r(0, i) * points3Dn(2, i);
				A(2 * i + 1, 1) = s(0, i) * points3Dn(2, i);

				//r22
				A(2 * i, 2) = r(1, i) * points3Dn(1, i);
				A(2 * i + 1, 2) = s(1, i) * points3Dn(1, i);

				//r23
				A(2 * i, 3) = r(1, i) * points3Dn(2, i);
				A(2 * i + 1, 3) = s(1, i) * points3Dn(2, i);

				//r31
				A(2 * i, 4) = r(2, i) * points3Dn(1, i);
				A(2 * i + 1, 4) = s(2, i) * points3Dn(1, i);

				//r32
				A(2 * i, 5) = r(2, i) * points3Dn(2, i);
				A(2 * i + 1, 5) = s(2, i) * points3Dn(2, i);

				//t1
				A(2 * i, 6) = r(0, i);
				A(2 * i + 1, 6) = s(0, i);

				//t2
				A(2 * i, 7) = r(1, i);
				A(2 * i + 1, 7) = s(1, i);

				//t3
				A(2 * i, 8) = r(2, i);
				A(2 * i + 1, 8) = s(2, i);

			}

		}
		else
		{
			// build full system.

			for (std::size_t i = 0; i < nrPts; i++)
			{

				//r11
				A(2 * i, 0) = r(0, i) * points3Dn(0, i);
				A(2 * i + 1, 0) = s(0, i) * points3Dn(0, i);

				//r12
				A(2 * i, 1) = r(0, i) * points3Dn(1, i);
				A(2 * i + 1, 1) = s(0, i) * points3Dn(1, i);

				//r13
				A(2 * i, 2) = r(0, i) * points3Dn(2, i);
				A(2 * i + 1, 2) = s(0, i) * points3Dn(2, i);

				//r21
				A(2 * i, 3) = r(1, i) * points3Dn(0, i);
				A(2 * i + 1, 3) = s(1, i) * points3Dn(0, i);

				//r22
				A(2 * i, 4) = r(1, i) * points3Dn(1, i);
				A(2 * i + 1, 4) = s(1, i) * points3Dn(1, i);

				//r23
				A(2 * i, 5) = r(1, i) * points3Dn(2, i);
				A(2 * i + 1, 5) = s(1, i) * points3Dn(2, i);

				//r31
				A(2 * i, 6) = r(2, i) * points3Dn(0, i);
				A(2 * i + 1, 6) = s(2, i) * points3Dn(0, i);

				//r32
				A(2 * i, 7) = r(2, i) * points3Dn(1, i);
				A(2 * i + 1, 7) = s(2, i) * points3Dn(1, i);

				//r33
				A(2 * i, 8) = r(2, i) * points3Dn(2, i);
				A(2 * i + 1, 8) = s(2, i) * points3Dn(2, i);

				//t1
				A(2 * i, 9) = r(0, i);
				A(2 * i + 1, 9) = s(0, i);

				//t2
				A(2 * i, 10) = r(1, i);
				A(2 * i + 1, 10) = s(1, i);

				//t3
				A(2 * i, 11) = r(2, i);
				A(2 * i + 1, 11) = s(2, i);

			}
		}

		//do least squares (A^t)PAx = 0

		Eigen::MatrixXd b = A.transpose() * A;

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(b, Eigen::ComputeFullV);

		Eigen::MatrixXd v1 = svd.matrixV();


		if (planar == 1)
		{
			Eigen::Vector3d tout1 = v1.rightCols<1>().tail<3>();

			Eigen::Matrix3d P = Eigen::Matrix3d::Zero();

			Eigen::VectorXd temp = v1.rightCols<1>().head<6>();

			P.col(1)(0) = temp(0);
			P.col(1)(1) = temp(2);
			P.col(1)(2) = temp(4);

			P.col(2)(0) = temp(1);
			P.col(2)(1) = temp(3);
			P.col(2)(2) = temp(5);

			double scalefact = std::sqrt(std::abs(P.col(1).norm() * P.col(2).norm()));

			P.col(0) = P.col(1).cross(P.col(2));

			P = P.transpose();

			//SVD to find the best rotation matrix in the Frobenius sense

			Eigen::JacobiSVD<Eigen::Matrix3d> svd(P, Eigen::ComputeFullU | Eigen::ComputeFullV);

			Eigen::Matrix3d U2 = svd.matrixU();
			Eigen::Matrix3d V2 = svd.matrixV();

			Eigen::Matrix3d R = U2 * V2.transpose();

			if (R.determinant() < 0)
			{
				R = -1.0 * R;
			}


			// rotate solution back(see paper)

			R = eigRot * R;

			//recover translation

			Eigen::Vector3d tout = tout1 / scalefact;

			R = -1.0 * R.transpose();

			Eigen::Matrix3d R1 = R;

			Eigen::Matrix3d R2;
			R2.leftCols<2>() = -R.leftCols<2>();
			R2.col(2) = R.col(2);

			std::vector<Eigen::Matrix4d> Ts(4);

			Ts[0].setZero();
			Ts[0].block<3, 3>(0, 0) = R1;
			Ts[0].col(3).head(3) = tout;
			Ts[0](3, 3) = 1;

			Ts[1].setZero();
			Ts[1].block<3, 3>(0, 0) = R1;
			Ts[1].col(3).head(3) = -tout;
			Ts[1](3, 3) = 1;

			Ts[2].setZero();
			Ts[2].block<3, 3>(0, 0) = R2;
			Ts[2].col(3).head(3) = tout;
			Ts[2](3, 3) = 1;

			Ts[3].setZero();
			Ts[3].block<3, 3>(0, 0) = R2;
			Ts[3].col(3).head(3) = -tout;
			Ts[3](3, 3) = 1;

			//find the best solution with 6 correspondences

			Eigen::Vector4d diff1;
			diff1.setZero();


			for (int te = 0; te < 6; te++)
			{
				for (int ba = 0; ba < 4; ba++)
				{
					Eigen::Vector4d temp = Eigen::Vector4d::Ones();
					temp.head<3>() = points3D.col(te);

					Eigen::Vector4d testres1 = Ts[ba] * temp;

					Eigen::Vector3d testres11 = testres1.head<3>() / (testres1.head<3>().norm());

					diff1(ba) = diff1(ba) + (1 - testres11.dot(v.col(te)));

				}// end for (ba).

			} // end for (te).

			int idx;

			diff1.minCoeff(&idx);

			T = Ts[idx]; //注意, T是4x4 矩阵.

		}
		else
		{

			Eigen::Vector3d tout1 = v1.rightCols<1>().tail<3>();

			Eigen::Matrix3d P;
			P.col(0) = v1.rightCols<1>().head<3>();
			P.col(1) = v1.rightCols<1>().segment<3>(3);//取向量中间的3个元素.
			P.col(2) = v1.rightCols<1>().segment<3>(6);

			double scalefact = std::pow(abs(P.col(0).norm() * P.col(1).norm() * P.col(2).norm()), 1.0 / 3.0);


			//SVD to find the best rotation matrix in the Frobenius sense

			Eigen::JacobiSVD<Eigen::Matrix3d> svd(P, Eigen::ComputeFullU | Eigen::ComputeFullV);

			Eigen::Matrix3d U2 = svd.matrixU();
			Eigen::Matrix3d V2 = svd.matrixV();

			Eigen::Matrix3d R = U2 * V2.transpose();

			if (R.determinant() < 0)
			{
				R = -R;
			}

			// recover translation

			Eigen::Vector3d tout = R * (tout1 / scalefact);

			Eigen::Matrix4d temp1;
			temp1.setZero();
			temp1.block<3, 3>(0, 0) = R;
			temp1.col(3).head<3>() = tout;
			temp1(3, 3) = 1;
			Eigen::Matrix4d T1 = temp1.inverse();


			Eigen::Matrix4d temp2;
			temp2.setZero();
			temp2.block<3, 3>(0, 0) = R;
			temp2.col(3).head<3>() = -tout;
			temp2(3, 3) = 1;
			Eigen::Matrix4d T2 = temp2.inverse();


			//find the best solution with 6 correspondences
			double diff1 = 0, diff2 = 0;

			for (int te = 0; te < 6; te++)
			{
				Eigen::Vector4d temp = Eigen::Vector4d::Ones();
				temp.head<3>() = points3D.col(te);
				Eigen::Vector4d test1 = T1 * temp;

				Eigen::Vector4d test2 = T2 * temp;


				Eigen::Vector3d testres1 = test1.head<3>() / (test1.head<3>().norm());
				Eigen::Vector3d testres2 = test2.head<3>() / (test2.head<3>().norm());

				diff1 = diff1 + (1 - testres1.dot(v.col(te)));
				diff2 = diff2 + (1 - testres2.dot(v.col(te)));

			}// end for (te).


			if (diff1 < diff2)
			{
				T = T1.block<3, 4>(0, 0);

			}
			else
			{
				T = T2.block<3, 4>(0, 0);

			}

		}

		OptimFlags optimFlags;
		optimFlags.epsP = 1e-6;
		optimFlags.epsF = 1e-6;
		optimFlags.maxit = 5;
		optimFlags.tau = 1e-4;

		Eigen::Matrix<double, 3, 4, Eigen::ColMajor> tempT = T, tempTT;

		int flag = OptimMLPnP_GN(tempT, points3D, r, s, Kll, optimFlags, tempTT, statistics);

		if (flag == -1)
		{
			std::cout << " Error happend !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;
		}

		//最终结果.
		T = tempTT;

		return 1;
	}

	/******************************************************/
	/* MLPnP solver with COV.
	/* points3D, 输入参数, 表示世界坐标系下的三维点.
	/* v, 输入参数, 图像点(齐次坐标表示)
	/* cov, 输入参数,表示协方差矩阵.
	/* T, 输出参数, 表示3x4矩阵.
	/* statistics, 输出参数,
	/*******************************************************/

	int CMLPnP::MLPnPSolver_COV(
		Eigen::Matrix3Xd& points3D, Eigen::Matrix3Xd& v,
		Eigen::MatrixXd& cov, Eigen::MatrixXd& T, Statistics& statistics)
	{



		std::size_t nrPts = points3D.cols();
		std::size_t vn = v.cols();

		if (nrPts == 0)
		{
			std::cout << "points3D is empty !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;
		}

		if (vn == 0)
		{
			std::cout << "v is empty !" << " at line: "
				<< __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;
		}

		if (nrPts != vn)
		{
			std::cout << "The size of points3D is not same as that of v !"
				<< " at line: " << __LINE__ << " in file: " << __FILE__ << std::endl;

			return -1;

		}



		//matrix of null space vectors r and s

		Eigen::Matrix3Xd r = Eigen::Matrix3Xd::Zero(3, nrPts);
		Eigen::Matrix3Xd s = Eigen::Matrix3Xd::Zero(3, nrPts);

		std::vector<Eigen::Matrix2d> cov_reduced(nrPts);

		// test planarity, only works well if the scene is really planar
		// quasi - planar won't work very well

		Eigen::Matrix3d S = points3D * points3D.transpose();

		Eigen::EigenSolver<Eigen::Matrix3d> solver(S);

		Eigen::MatrixXd eigRot = solver.pseudoEigenvectors();

		int planar = 0;

		//create full design matrix

		Eigen::MatrixXd points3Dn;

		Eigen::MatrixXd A;

		//注意矩阵秩的求法.
		Eigen::FullPivLU <Eigen::Matrix3d> lu_decomp(S);

		if (lu_decomp.rank() == 2)
		{
			planar = 1;

			Eigen::Matrix2Xd points3D1 = eigRot.transpose() * points3D;

			points3Dn = Eigen::MatrixXd::Ones(3, nrPts);

			points3Dn.topRows<2>() = points3D1;

			//points3Dn 的第三行元素全为1.

			// create reduced design matrix

			A = Eigen::MatrixXd::Zero(2 * nrPts, 9);


		}
		else
		{

			points3Dn = Eigen::MatrixXd::Ones(4, nrPts);

			points3Dn.topRows<3>() = points3D;

			A = Eigen::MatrixXd::Zero(2 * nrPts, 12);

		}


		// compute null spaces of bearing vector v : null(v')

		//Eigen::MatrixXd null_2d = Eigen::MatrixXd::Zero(nrPts, 2);

		Eigen::MatrixXd null_2d;

		for (int i = 0; i < nrPts; i++)
		{

			//注意向量核空间的求解方法.
			// nullspace of right vector.
			Eigen::JacobiSVD<Eigen::MatrixXd, Eigen::HouseholderQRPreconditioner>
				svd_solver(v.col(i).transpose(), Eigen::ComputeFullV);

			//null_2d = svd_solver.matrixV().block(0, 1, 3, 2);
			null_2d = svd_solver.matrixV().rightCols<2>();

			r.col(i) = null_2d.col(0);
			s.col(i) = null_2d.col(1);

			// if use_cov
			Eigen::Matrix3d tmp;
			tmp(0, 0) = cov.col(0)(0);
			tmp(1, 0) = cov.col(0)(1);
			tmp(2, 0) = cov.col(0)(2);

			tmp(0, 1) = cov.col(0)(3);
			tmp(1, 1) = cov.col(0)(4);
			tmp(2, 1) = cov.col(0)(5);

			tmp(0, 2) = cov.col(0)(6);
			tmp(1, 2) = cov.col(0)(7);
			cov(2, 2) = cov.col(0)(8);

			cov_reduced[i] = (null_2d.transpose() * tmp * null_2d).inverse();

		}

		// stochastic model

		Eigen::MatrixXd Kll = Eigen::MatrixXd::Identity(2 * nrPts, 2 * nrPts);

		if (planar == 1)
		{
			// build reduces system.


			for (std::size_t i = 0; i < nrPts; i++)
			{
				// if use_cov
				Kll.block(2 * i, 2 * i + 1, 2, 2) = cov_reduced[i];


				//r12.			
				A(2 * i, 0) = r(0, i) * points3Dn(1, i);
				A(2 * i + 1, 0) = s(0, i) * points3Dn(1, i);

				//r13
				A(2 * i, 1) = r(0, i) * points3Dn(2, i);
				A(2 * i + 1, 1) = s(0, i) * points3Dn(2, i);

				//r22
				A(2 * i, 2) = r(1, i) * points3Dn(1, i);
				A(2 * i + 1, 2) = s(1, i) * points3Dn(1, i);

				//r23
				A(2 * i, 3) = r(1, i) * points3Dn(2, i);
				A(2 * i + 1, 3) = s(1, i) * points3Dn(2, i);

				//r31
				A(2 * i, 4) = r(2, i) * points3Dn(1, i);
				A(2 * i + 1, 4) = s(2, i) * points3Dn(1, i);

				//r32
				A(2 * i, 5) = r(2, i) * points3Dn(2, i);
				A(2 * i + 1, 5) = s(2, i) * points3Dn(2, i);

				//t1
				A(2 * i, 6) = r(0, i);
				A(2 * i + 1, 6) = s(0, i);

				//t2
				A(2 * i, 7) = r(1, i);
				A(2 * i + 1, 7) = s(1, i);

				//t3
				A(2 * i, 8) = r(2, i);
				A(2 * i + 1, 8) = s(2, i);

			}

		}
		else
		{
			// build full system.

			for (std::size_t i = 0; i < nrPts; i++)
			{
				// if use_cov
				Kll.block(2 * i, 2 * i + 1, 2, 2) = cov_reduced[i];

				//r11
				A(2 * i, 0) = r(0, i) * points3Dn(0, i);
				A(2 * i + 1, 0) = s(0, i) * points3Dn(0, i);

				//r12
				A(2 * i, 1) = r(0, i) * points3Dn(1, i);
				A(2 * i + 1, 1) = s(0, i) * points3Dn(1, i);

				//r13
				A(2 * i, 2) = r(0, i) * points3Dn(2, i);
				A(2 * i + 1, 2) = s(0, i) * points3Dn(2, i);

				//r21
				A(2 * i, 3) = r(1, i) * points3Dn(0, i);
				A(2 * i + 1, 3) = s(1, i) * points3Dn(0, i);

				//r22
				A(2 * i, 4) = r(1, i) * points3Dn(1, i);
				A(2 * i + 1, 4) = s(1, i) * points3Dn(1, i);

				//r23
				A(2 * i, 5) = r(1, i) * points3Dn(2, i);
				A(2 * i + 1, 5) = s(1, i) * points3Dn(2, i);

				//r31
				A(2 * i, 6) = r(2, i) * points3Dn(0, i);
				A(2 * i + 1, 6) = s(2, i) * points3Dn(0, i);

				//r32
				A(2 * i, 7) = r(2, i) * points3Dn(1, i);
				A(2 * i + 1, 7) = s(2, i) * points3Dn(1, i);

				//r33
				A(2 * i, 8) = r(2, i) * points3Dn(2, i);
				A(2 * i + 1, 8) = s(2, i) * points3Dn(2, i);

				//t1
				A(2 * i, 9) = r(0, i);
				A(2 * i + 1, 9) = s(0, i);

				//t2
				A(2 * i, 10) = r(1, i);
				A(2 * i + 1, 10) = s(1, i);

				//t3
				A(2 * i, 11) = r(2, i);
				A(2 * i + 1, 11) = s(2, i);

			}
		}

		//do least squares (A^t)PAx = 0

		Eigen::MatrixXd b = A.transpose() * A;

		Eigen::JacobiSVD<Eigen::MatrixXd> svd(b, Eigen::ComputeFullV);

		Eigen::MatrixXd v1 = svd.matrixV();


		if (planar == 1)
		{
			Eigen::Vector3d tout1 = v1.rightCols<1>().tail<3>();

			Eigen::Matrix3d P = Eigen::Matrix3d::Zero();

			Eigen::VectorXd temp = v1.rightCols<1>().head<6>();

			P.col(1)(0) = temp(0);
			P.col(1)(1) = temp(2);
			P.col(1)(2) = temp(4);

			P.col(2)(0) = temp(1);
			P.col(2)(1) = temp(3);
			P.col(2)(2) = temp(5);

			double scalefact = std::sqrt(std::abs(P.col(1).norm() * P.col(2).norm()));

			P.col(0) = P.col(1).cross(P.col(2));

			P = P.transpose();

			//SVD to find the best rotation matrix in the Frobenius sense

			Eigen::JacobiSVD<Eigen::Matrix3d> svd(P, Eigen::ComputeFullU | Eigen::ComputeFullV);

			Eigen::Matrix3d U2 = svd.matrixU();
			Eigen::Matrix3d V2 = svd.matrixV();

			Eigen::Matrix3d R = U2 * V2.transpose();

			if (R.determinant() < 0)
			{
				R = -1.0 * R;
			}


			// rotate solution back(see paper)

			R = eigRot * R;

			//recover translation

			Eigen::Vector3d tout = tout1 / scalefact;

			R = -1.0 * R.transpose();

			Eigen::Matrix3d R1 = R;

			Eigen::Matrix3d R2;
			R2.leftCols<2>() = -R.leftCols<2>();
			R2.col(2) = R.col(2);

			std::vector<Eigen::Matrix4d> Ts(4);

			Ts[0].setZero();
			Ts[0].block<3, 3>(0, 0) = R1;
			Ts[0].col(3).head(3) = tout;
			Ts[0](3, 3) = 1;

			Ts[1].setZero();
			Ts[1].block<3, 3>(0, 0) = R1;
			Ts[1].col(3).head(3) = -tout;
			Ts[1](3, 3) = 1;

			Ts[2].setZero();
			Ts[2].block<3, 3>(0, 0) = R2;
			Ts[2].col(3).head(3) = tout;
			Ts[2](3, 3) = 1;

			Ts[3].setZero();
			Ts[3].block<3, 3>(0, 0) = R2;
			Ts[3].col(3).head(3) = -tout;
			Ts[3](3, 3) = 1;

			//find the best solution with 6 correspondences

			Eigen::Vector4d diff1;
			diff1.setZero();


			for (int te = 0; te < 6; te++)
			{
				for (int ba = 0; ba < 4; ba++)
				{
					Eigen::Vector4d temp = Eigen::Vector4d::Ones();
					temp.head<3>() = points3D.col(te);

					Eigen::Vector4d testres1 = Ts[ba] * temp;

					Eigen::Vector3d testres11 = testres1.head<3>() / (testres1.head<3>().norm());

					diff1(ba) = diff1(ba) + (1 - testres11.dot(v.col(te)));

				}// end for (ba).

			} // end for (te).

			int idx;

			diff1.minCoeff(&idx);

			T = Ts[idx]; //注意, T是4x4 矩阵.

		}
		else
		{

			Eigen::Vector3d tout1 = v1.rightCols<1>().tail<3>();

			Eigen::Matrix3d P;
			P.col(0) = v1.rightCols<1>().head<3>();
			P.col(1) = v1.rightCols<1>().segment<3>(3);//取向量中间的3个元素.
			P.col(2) = v1.rightCols<1>().segment<3>(6);

			double scalefact = std::pow(abs(P.col(0).norm() * P.col(1).norm() * P.col(2).norm()), 1.0 / 3.0);


			//SVD to find the best rotation matrix in the Frobenius sense

			Eigen::JacobiSVD<Eigen::Matrix3d> svd(P, Eigen::ComputeFullU | Eigen::ComputeFullV);

			Eigen::Matrix3d U2 = svd.matrixU();
			Eigen::Matrix3d V2 = svd.matrixV();

			Eigen::Matrix3d R = U2 * V2.transpose();

			if (R.determinant() < 0)
			{
				R = -R;
			}

			// recover translation

			Eigen::Vector3d tout = R * (tout1 / scalefact);

			Eigen::Matrix4d temp1;
			temp1.setZero();
			temp1.block<3, 3>(0, 0) = R;
			temp1.col(3).head<3>() = tout;
			temp1(3, 3) = 1;
			Eigen::Matrix4d T1 = temp1.inverse();


			Eigen::Matrix4d temp2;
			temp2.setZero();
			temp2.block<3, 3>(0, 0) = R;
			temp2.col(3).head<3>() = -tout;
			temp2(3, 3) = 1;
			Eigen::Matrix4d T2 = temp2.inverse();


			//find the best solution with 6 correspondences
			double diff1 = 0, diff2 = 0;

			for (int te = 0; te < 6; te++)
			{
				Eigen::Vector4d temp = Eigen::Vector4d::Ones();
				temp.head<3>() = points3D.col(te);
				Eigen::Vector4d test1 = T1 * temp;

				Eigen::Vector4d test2 = T2 * temp;


				Eigen::Vector3d testres1 = test1.head<3>() / (test1.head<3>().norm());
				Eigen::Vector3d testres2 = test2.head<3>() / (test2.head<3>().norm());

				diff1 = diff1 + (1 - testres1.dot(v.col(te)));
				diff2 = diff2 + (1 - testres2.dot(v.col(te)));

			}// end for (te).


			if (diff1 < diff2)
			{
				T = T1.block<3, 4>(0, 0);

			}
			else
			{
				T = T2.block<3, 4>(0, 0);

			}

		}

		OptimFlags optimFlags;
		optimFlags.epsP = 1e-6;
		optimFlags.epsF = 1e-6;
		optimFlags.maxit = 5;
		optimFlags.tau = 1e-4;

		Eigen::Matrix<double, 3, 4, Eigen::ColMajor> tempT = T, tempTT;

		int flag = OptimMLPnP_GN(tempT, points3D, r, s, Kll, optimFlags, tempTT, statistics);

		if (flag == -1)
		{
			std::cout << " Error happend !" << " at line: "
				<< __LINE__ << ", in file: " << __FILE__ << std::endl;

			return -1;
		}

		//最终结果.
		T = tempTT;


		return 1;

	}

}//namespace cvg;
