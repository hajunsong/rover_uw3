#pragma once

#include <Eigen/Core>
#include <Eigen/Dense>
#include <iostream>
#include <vector>

#include "fileio.h"

class Body {
public:
	Body(){u_vec << 0, 0, 1;};
	~Body(){};

	int id;

	// intput data
	double qi, dqi, ddqi, mi;
	Eigen::Vector3d ri, rhoip, sijp;
	Eigen::Matrix3d Jip, Cii, Cij;
	Eigen::Vector3d u_vec;

	// position anlaysis
	Eigen::Matrix3d Aijpp, Ai;
	Eigen::Vector3d sij, rhoi, ric;

	// vleocity analysis
	Eigen::Vector3d Hi, wi, dri, dric, dHi;
	Eigen::Matrix3d wit, rit, drit;
	Eigen::VectorXd Bi, Di, Yih;

	// mass force analysis
	Eigen::Vector3d fic, tic;
	Eigen::Matrix3d Jic, Ai_Cii, rict, drict;
	Eigen::MatrixXd Mih, Ki;
	Eigen::VectorXd Qih, Li;

	// acceleration analysis
	Eigen::VectorXd dYih, Ri, dYib;
	Eigen::MatrixXd Ti, dTi;
	Eigen::Vector3d ddri, dwi, ddric;
	Eigen::Matrix3d dwit;

	// inverse dynamics
	Eigen::VectorXd Rjih, Rjib;
	Eigen::Vector3d fji, nji, fjip, njip;
};

class RobotArm {
public:
	RobotArm(unsigned char NumBody, unsigned char DOF);
	~RobotArm();

	void run();

	unsigned char num_body, dof;

private:
	std::vector<Body> body;

	double t_current, t_end, g, h;
	FILE* fp;

	Eigen::Matrix3d tilde(Eigen::Vector3d x);

	void read_base(Body* body);
	void read_body1(Body* body);
	void read_body2(Body* body);
	void read_body3(Body* body);
	void read_body4(Body* body);
	void read_body5(Body* body);
	void read_body6(Body* body);

	void position_analysis(Body *jbody, Body *ibody);
	void velocity_analysis(Body *jbody, Body *ibody);
	void mass_force_analysis(Body *ibody);
	void acceleration_analysis(Body *jbody, Body *ibody);
	void inverse_dynamics(Body *ibody);

	void data_save();
};
