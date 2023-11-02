#include "robotarm.h"

RobotArm::RobotArm(unsigned char NumBody, unsigned char DOF)
{
	num_body = NumBody;
	dof = DOF;
	body.resize(num_body + 1);

	std::cout << std::setprecision(20) << std::fixed;
}

RobotArm::~RobotArm()
{
	body.clear();
}

void RobotArm::run()
{
	std::string q_file_name = "../data/recurdyn_q.csv";
	std::vector<double> q_data;
	uint q_data_size[2] = {0, 0};
	uint q_row, q_col;
	load_data(q_file_name, &q_data, ",", q_data_size);
	q_row = q_data_size[0];
	q_col = q_data_size[1];

	std::string f_file_name = "../data/recurdyn_force.csv";
	std::vector<double> f_data;
	uint f_data_size[2] = {0, 0};
	uint f_row, f_col;
	load_data(f_file_name, &f_data, ",", f_data_size);
	f_row = f_data_size[0];
	f_col = f_data_size[1];

	sensor_force = Eigen::VectorXd(6);

	h = q_data[2*q_col + 0] - q_data[1*q_col + 0];
	g = -9.80665;

	read_base(&body[0]);
	read_body1(&body[1]);
	read_body2(&body[2]);
	read_body3(&body[3]);
	read_body4(&body[4]);
	read_body5(&body[5]);
	read_body6(&body[6]);

	fp = fopen("../data/cpp_torque.txt", "w+");

	for(uint indx = 3; indx < q_row; indx++){
		t_current = q_data[indx*q_col + 0];

		for(uint i = 1; i <= num_body; i++){
			body[i].qi = q_data[indx*q_col + i];
			body[i].dqi = (q_data[indx*q_col + i] - q_data[(indx - 1)*q_col + i])/h;
			body[i].ddqi = (((q_data[indx*q_col + i] - q_data[(indx - 1)*q_col + i])/h) - ((q_data[(indx - 1)*q_col + i] - q_data[(indx - 2)*q_col + i])/h))/h;
		}

		for(uint i = 0; i < 6; i++){
			sensor_force(i) = f_data[indx*f_col + i + 1];
		}

		for(uint i = 1; i <= num_body; i++){
			position_analysis(&body[i], &body[i - 1]);
			velocity_analysis(&body[i], &body[i - 1]);
		}

		for(uint i = num_body; i > 0; i--){
			mass_force_analysis(&body[i]);
		}

		for(uint i = 1; i <= num_body; i++){
			acceleration_analysis(&body[i], &body[i - 1]);
		}

		for(uint i = num_body; i > 0; i--){
			inverse_dynamics(&body[i]);
		}

		data_save();

		std::cout << "t_current : " << t_current << std::endl;
	}

	fclose(fp);
}

void RobotArm::position_analysis(Body *jbody, Body *ibody)
{
	jbody->Aijpp(0,0) = cos(jbody->qi);	jbody->Aijpp(0,1) = -sin(jbody->qi);jbody->Aijpp(0,2) = 0;
	jbody->Aijpp(1,0) = sin(jbody->qi);	jbody->Aijpp(1,1) = cos(jbody->qi);	jbody->Aijpp(1,2) = 0;
	jbody->Aijpp(2,0) = 0;				jbody->Aijpp(2,1) = 0;				jbody->Aijpp(2,2) = 1;

	jbody->Ai = ibody->Ai*ibody->Cij*jbody->Aijpp;
	jbody->sij = jbody->Ai*jbody->sijp;
	jbody->ri = ibody->ri + ibody->sij;
	jbody->rhoi = jbody->Ai*jbody->rhoip;
	jbody->ric = jbody->ri + jbody->rhoi;
}

void RobotArm::velocity_analysis(Body *jbody, Body *ibody)
{
	jbody->Hi = ibody->Ai*ibody->Cij*jbody->u_vec;
	jbody->wi = ibody->wi + jbody->Hi*jbody->dqi;
	jbody->wit = tilde(jbody->wi);
	jbody->dri = ibody->dri + ibody->wit*ibody->sij;
	jbody->rit = tilde(jbody->ri);

	jbody->Bi = Eigen::VectorXd(6);
	jbody->Bi << jbody->rit*jbody->Hi
					 , jbody->Hi;

	jbody->drit= tilde(jbody->dri);
	jbody->dric= jbody->dri+ jbody->wit*jbody->rhoi;

	jbody->dHi= ibody->wit*jbody->Hi;

	jbody->Di = Eigen::VectorXd(6);
	jbody->Di << jbody->drit*jbody->Hi + jbody->rit*jbody->dHi
					 , jbody->dHi;
	jbody->Di *= jbody->dqi;

	jbody->Yih = ibody->Yih + jbody->Bi*jbody->dqi;
}

void RobotArm::mass_force_analysis(Body *ibody)
{
	ibody->Ai_Cii = ibody->Ai*ibody->Cii;
	ibody->Jic = ibody->Ai_Cii*ibody->Jip*ibody->Ai_Cii.transpose();

	ibody->rict = tilde(ibody->ric);
	ibody->drict = tilde(ibody->dric);

	ibody->fic << 0, 0, ibody->mi*g;
	ibody->tic << 0, 0, 0;

	if(ibody->id == num_body){
		ibody->fic += sensor_force.segment(0, 3);
		ibody->tic += sensor_force.segment(3, 3)*0.001;
	}

	ibody->Mih = Eigen::MatrixXd(6,6);
	ibody->Mih.block(0, 0, 3, 3) = ibody->mi*Eigen::Matrix3d::Identity();
	ibody->Mih.block(0, 3, 3, 3) = -ibody->mi*ibody->rict;
	ibody->Mih.block(3, 0, 3, 3) = -ibody->Mih.block(0, 3, 3, 3);
	ibody->Mih.block(3, 3, 3, 3) = ibody->Jic - ibody->mi*ibody->rict*ibody->rict;

	ibody->Qih = Eigen::VectorXd(6);
	ibody->Qih << ibody->fic + ibody->mi*ibody->drict*ibody->wi
					, ibody->tic + ibody->rict*ibody->fic + ibody->mi*ibody->rict*ibody->drict*ibody->wi - ibody->wit*ibody->Jic*ibody->wi;

	ibody->Ki = Eigen::MatrixXd(6, 6);
	ibody->Ki = ibody->Mih;

	ibody->Li = Eigen::VectorXd(6);
	ibody->Li = ibody->Qih;

	if(ibody->id != num_body){
		Body *jbody = &body[ibody->id + 1];
		ibody->Ki += jbody->Ki;
		ibody->Li += jbody->Li - jbody->Ki*jbody->Di;
	}
}

void RobotArm::acceleration_analysis(Body *jbody, Body *ibody)
{
	jbody->dYih = Eigen::VectorXd(6);
	jbody->dYih = ibody->dYih + jbody->Bi*jbody->ddqi + jbody->Di;

	jbody->Ti = Eigen::MatrixXd(6, 6);
	jbody->Ti.block(0, 0, 3, 3) = Eigen::MatrixXd::Identity(3, 3);
	jbody->Ti.block(0, 3, 3, 3) = Eigen::MatrixXd::Zero(3, 3);
	jbody->Ti.block(3, 0, 3, 3) = -jbody->rit;
	jbody->Ti.block(3, 3, 3, 3) = jbody->Ti.block(0, 0, 3, 3);

	jbody->dTi = Eigen::MatrixXd::Zero(6,6);
	jbody->dTi.block(0, 3, 3, 3) = -jbody->drit;

	jbody->Ri = Eigen::VectorXd::Zero(6);
	jbody->Ri.segment(0, 3) = jbody->drit*jbody->wi;

	jbody->dYib = Eigen::VectorXd(6);
	jbody->dYib = jbody->dTi*jbody->Yih + jbody->Ti*jbody->dYih;

	jbody->ddri = jbody->dYib.segment(0, 3);
	jbody->dwi = jbody->dYib.segment(3, 3);

	jbody->dwit = tilde(jbody->dwi);

	jbody->ddric = jbody->ddri + jbody->dwit*jbody->rhoi + jbody->wit*jbody->wit*jbody->rhoi;
}

void RobotArm::inverse_dynamics(Body *ibody)
{
	ibody->Rjih = Eigen::VectorXd(6);
	ibody->Rjih = -ibody->Li + ibody->Ki*ibody->dYih;
	for(uint i = ibody->id + 1; i <= num_body; i++){
		ibody->Rjih += body[i].Ki*body[i].Bi*body[i].ddqi;
	}

	ibody->Rjib = ibody->Ti*ibody->Rjih;

	ibody->fji = ibody->Rjib.segment(0, 3);
	ibody->nji = ibody->Rjib.segment(3, 3);

	ibody->fjip = ibody->Ai.transpose()*ibody->fji;
	ibody->njip = ibody->Ai.transpose()*ibody->nji;

	// std::cout << "id : " << ibody->id << std::endl;
	// std::cout << ibody->Ti << std::endl << std::endl;
}

void RobotArm::data_save()
{
	fprintf(fp, "%3.5f ", t_current);
	for(uint i = 1; i <= num_body; i++){
		fprintf(fp, "%3.7f ", body[i].njip[2]);
	}
	fprintf(fp, "\n");
}

Eigen::Matrix3d RobotArm::tilde(Eigen::Vector3d x)
{
	Eigen::Matrix3d x_tilde;

	x_tilde(0, 0) = 0.0;	x_tilde(0, 1) = -x(2);	x_tilde(0, 2) = x(1);
	x_tilde(1, 0) = x(2);	x_tilde(1, 1) = 0.0;	x_tilde(1, 2) = -x(0);
	x_tilde(2, 0) = -x(1);	x_tilde(2, 1) = x(0);	x_tilde(2, 2) = 0.0;

	return x_tilde;
}

void RobotArm::read_base(Body *body)
{
	body->id = 0;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 734.81362;

	body->Jip(0,0) = 100.95495; body->Jip(0,1) = 0;         body->Jip(0,2) = 0;
	body->Jip(1,0) = 0;			body->Jip(1,1) = 133.20808;	body->Jip(1,2) = 0;
	body->Jip(2,0) = 0;			body->Jip(2,1) = 0;         body->Jip(2,2) = 140.46847;

	body->rhoip(0) = 0.3;
	body->rhoip(1) = 0;
	body->rhoip(2) = 0.261986;

	body->Cii = Eigen::Matrix3d::Identity();
	body->Cij = Eigen::Matrix3d::Identity();

	body->sijp(0) = 0.3;
	body->sijp(1) = 0;
	body->sijp(2) = 0.4508509;

	body->ri(0) = 0;
	body->ri(1) = 0;
	body->ri(2) = 0.5976491;

	body->dri = Eigen::Vector3d::Zero();
	body->wi = Eigen::Vector3d::Zero();

	body->Ai = Eigen::Matrix3d::Identity();
	body->sij = body->Ai*body->sijp;
	body->rhoi =  body->Ai*body->rhoip;
	body->ric = body->ri + body->rhoi;

	body->wit = tilde(body->wi);
	body->Yih = Eigen::VectorXd(6);
	body->drit= tilde(body->dri);
	body->Yih << body->dri + body->drit*body->wi
				 , body->wi;
	body->dYih = Eigen::VectorXd::Zero(6);
}

void RobotArm::read_body1(Body *body)
{
	body->id = 1;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 38.6603;

	body->Jip(0,0) = 0.25617402; body->Jip(0,1) = 0;          body->Jip(0,2) = 0;
	body->Jip(1,0) = 0;          body->Jip(1,1) = 0.41035289; body->Jip(1,2) = 0;
	body->Jip(2,0) = 0;          body->Jip(2,1) = 0;          body->Jip(2,2) = 0.32227293;

	body->sijp(0) = 0.171;
	body->sijp(1) = 0;
	body->sijp(2) = 0;

	body->Cij(0,0) = 1; body->Cij(0,1) = 0;  body->Cij(0,2) = 0;
	body->Cij(1,0) = 0; body->Cij(1,1) = 0;  body->Cij(1,2) = 1;
	body->Cij(2,0) = 0; body->Cij(2,1) = -1; body->Cij(2,2) = 0;

	body->rhoip(0) = 0.0560604;
	body->rhoip(1) = -5.23142e-04;
	body->rhoip(2) = -0.0719175;

	body->Cii(0,0) = 1; body->Cii(0,1) = 0; body->Cii(0,2) = 0;
	body->Cii(1,0) = 0; body->Cii(1,1) = 1; body->Cii(1,2) = 0;
	body->Cii(2,0) = 0; body->Cii(2,1) = 0; body->Cii(2,2) = 1;
}

void RobotArm::read_body2(Body *body)
{
	body->id = 2;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 110.76523;

	body->Jip(0,0) = 0.94029053; body->Jip(0,1) = 0;         body->Jip(0,2) = 0;
	body->Jip(1,0) = 0;          body->Jip(1,1) = 8.4072768; body->Jip(1,2) = 0;
	body->Jip(2,0) = 0;          body->Jip(2,1) = 0;         body->Jip(2,2) = 8.588407;

	body->sijp(0) = 0.9213577;
	body->sijp(1) = 4.827e-07;
	body->sijp(2) = 1.46e-11;

	body->Cij(0,0) = 1; body->Cij(0,1) = 0;  body->Cij(0,2) = 0;
	body->Cij(1,0) = 0; body->Cij(1,1) = 1;  body->Cij(1,2) = 0;
	body->Cij(2,0) = 0; body->Cij(2,1) = 0;  body->Cij(2,2) = 1;

	body->rhoip(0) = 0.473039;
	body->rhoip(1) = -0.0162355;
	body->rhoip(2) = -0.0220838;

	body->Cii(0,0) = 0.798341020580324; body->Cii(0,1) = 2.21202477524923e-06;	body->Cii(0,2) = 0.602205625059974;
	body->Cii(1,0) = 0.602205625064037; body->Cii(1,1) = -2.93247031100656e-06; body->Cii(1,2) = -0.798341020574938;
	body->Cii(2,0) = 0;					body->Cii(2,1) = 0.999999999993254;		body->Cii(2,2) = -3.67320510334657e-06;
}

void RobotArm::read_body3(Body *body)
{
	body->id = 3;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 72.734409;

	body->Jip(0,0) = 2.2925008; body->Jip(0,1) = 0;         body->Jip(0,2) = 0;
	body->Jip(1,0) = 0;			body->Jip(1,1) = 2.2332289; body->Jip(1,2) = 0;
	body->Jip(2,0) = 0;			body->Jip(2,1) = 0;         body->Jip(2,2) = 0.40505243;

	body->sijp(0) = 0.535943381;
	body->sijp(1) = -4.98e-08;
	body->sijp(2) = -5.03e-12;

	body->Cij(0,0) = 1; body->Cij(0,1) = 0;  body->Cij(0,2) = 0;
	body->Cij(1,0) = 0; body->Cij(1,1) = 1;  body->Cij(1,2) = 0;
	body->Cij(2,0) = 0; body->Cij(2,1) = 0;  body->Cij(2,2) = 1;

	body->rhoip(0) = 0.535942;
	body->rhoip(1) = 0.00111722;
	body->rhoip(2) = -1.46e-11;

	body->Cii(0,0) = 0.769617146910743;	 body->Cii(0,1) = -2.34536216060518e-06; body->Cii(0,2) = -0.638505635977841;
	body->Cii(1,0) = -0.638505635982148; body->Cii(1,1) = -2.82696163165557e-06; body->Cii(1,2) = -0.769617146905551;
	body->Cii(2,0) = 0;					 body->Cii(2,1) = 0.999999999993254;	 body->Cii(2,2) = -3.67320510334657e-06;
}

void RobotArm::read_body4(Body *body)
{
	body->id = 4;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 29.19895;

	body->Jip(0,0) = 0.1627909; body->Jip(0,1) = 0;          body->Jip(0,2) = 0;
	body->Jip(1,0) = 0;			body->Jip(1,1) = 0.23597346; body->Jip(1,2) = 0;
	body->Jip(2,0) = 0;			body->Jip(2,1) = 0;          body->Jip(2,2) = 0.25574398;

	body->sijp(0) = 0.146;
	body->sijp(1) = 7.17e-12;
	body->sijp(2) = 3.28e-11;

	body->Cij(0,0) = 0; body->Cij(0,1) = 1;  body->Cij(0,2) = 0;
	body->Cij(1,0) = 0; body->Cij(1,1) = 0;  body->Cij(1,2) = 1;
	body->Cij(2,0) = 1; body->Cij(2,1) = 0;  body->Cij(2,2) = 0;

	body->rhoip(0) = 0.0715249;
	body->rhoip(1) = 0.0149249;
	body->rhoip(2) = -0.0145652;

	body->Cii(0,0) = 1;					body->Cii(0,1) = 0;						body->Cii(0,2) = 0;
	body->Cii(1,0) = 0;					body->Cii(1,1) = -3.67320510334657e-06; body->Cii(1,2) = -0.999999999993254;
	body->Cii(2,0) = 0;					body->Cii(2,1) = 0.999999999993254;		body->Cii(2,2) = -3.67320510334657e-06;
}

void RobotArm::read_body5(Body *body)
{
	body->id = 5;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 20.98815;

	body->Jip(0,0) = 0.093006759; body->Jip(0,1) = 0;          body->Jip(0,2) = 0;
	body->Jip(1,0) = 0;			  body->Jip(1,1) = 0.19180177; body->Jip(1,2) = 0;
	body->Jip(2,0) = 0;			  body->Jip(2,1) = 0;          body->Jip(2,2) = 0.14773127;

	body->sijp(0) = 5.53e-13;
	body->sijp(1) = 0.238;
	body->sijp(2) = -1.3e-11;

	body->Cij(0,0) = 0;  body->Cij(0,1) = -1; body->Cij(0,2) = 0;
	body->Cij(1,0) = 0;	 body->Cij(1,1) = 0;  body->Cij(1,2) = 1;
	body->Cij(2,0) = -1; body->Cij(2,1) = 0;  body->Cij(2,2) = 0;

	body->rhoip(0) = 0.000573087;
	body->rhoip(1) = 0.0643968;
	body->rhoip(2) = -0.00891077;

	body->Cii(0,0) = -3.67320510334657e-06;	body->Cii(0,1) = 0.999999999989733;		body->Cii(0,2) = 2.65358979333483e-06;
	body->Cii(1,0) = 0.999999999993254;		body->Cii(1,1) = 3.67320510333364e-06;	body->Cii(1,2) = 9.74717957113163e-12;
	body->Cii(2,0) = 0;						body->Cii(2,1) = 2.65358979335273e-06;	body->Cii(2,2) = -0.999999999996479;
}

void RobotArm::read_body6(Body *body)
{
	body->id = 6;

	body->qi = 0;
	body->dqi= 0;
	body->ddqi= 0;

	body->mi = 70.745182;

	body->Jip(0,0) = 0.687119280924208;		body->Jip(0,1) = -0.055977537779502;	body->Jip(0,2) = -1.07220027356336e-05;
	body->Jip(1,0) = -0.055977537779502;	body->Jip(1,1) = 1.18118135417106;		body->Jip(1,2) = -6.68198238753321e-06;
	body->Jip(2,0) = -1.07220027356336e-05;	body->Jip(2,1) = -6.68198238753321e-06; body->Jip(2,2) = 1.05307348253552;

	body->sijp(0) = -0.000858527;
	body->sijp(1) = -0.00176628;
	body->sijp(2) = 0.199286;

	body->Cij(0,0) = 0;	body->Cij(0,1) = 1;	body->Cij(0,2) = 0;
	body->Cij(1,0) = 0;	body->Cij(1,1) = 0;	body->Cij(1,2) = 1;
	body->Cij(2,0) = 1;	body->Cij(2,1) = 0;	body->Cij(2,2) = 0;

	body->rhoip(0) = -0.000865807;
	body->rhoip(1) = 0.0362375;
	body->rhoip(2) = 0.138304;

	body->Cii(0,0) = 3.67321859573274e-06;	body->Cii(0,1) = 3.67319161088606e-06;	body->Cii(0,2) = 0.999999999986507;
	body->Cii(1,0) = -3.67319161088606e-06;	body->Cii(1,1) = -0.999999999986507;	body->Cii(1,2) = 3.67320510332179e-06;
	body->Cii(2,0) = 0.999999999986507;		body->Cii(2,1) = -3.67320510332179e-06;	body->Cii(2,2) = -3.67320510334657e-06;
}

