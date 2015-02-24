#include "header.h"

void DisplayInfo(void){
	std::cout << "Rank of this process is " << PS::Comm::getRank() << std::endl;
	std::cout << "# of proc is " << PS::Comm::getNumberOfProc() << std::endl;
	return ;
}
