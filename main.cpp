//#include <igl/readOBJ.h>
#include <string>
#include <vector>
#include "Scripts.h"
#include <Fragment.h>
#include <Segment.h>


int main(int argc, char* argv[])
{
	std::string current_exec_name = argv[0]; // Name of the current exec program
	std::vector<std::string> all_args;

	if (argc > 1) {
		all_args.assign(argv + 1, argv + argc);
	}

	std::string fragmentPath = all_args[0];
	ObjFragment fragment = ObjFragment(fragmentPath);
	std::cout << "Loading data" << std::endl;
	fragment.load();

	std::string outFileName = all_args[1];
	std::string isSave = all_args[2];

	//segment_intact_surface(fragment,outFileName);
	//segment_opposite_surface(fragment);
	segment_sidewalls_surface(fragment);
}
