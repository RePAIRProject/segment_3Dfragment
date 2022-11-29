//#include <igl/readOBJ.h>
#include <string>
#include <vector>
#include "Scripts.h"



int main(int argc, char* argv[])
{
	std::string current_exec_name = argv[0]; // Name of the current exec program
	std::vector<std::string> all_args;

	if (argc > 1) {
		all_args.assign(argv + 1, argv + argc);
	}


	//std::string fragmentName = "RPf_00154"; //"cube";
	//std::string fragmentPath = "..\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";
	//std::string outputPath = "..\\fragments\\" + fragmentName + "\\" ;

	segment_intact_surface(all_args);
}
