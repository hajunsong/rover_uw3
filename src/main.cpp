#include "robotarm.h"

int main(){
	RobotArm robotArm(6, 6);

	std::cout << "start program" << std::endl;
	robotArm.run();

	return 0;
}
