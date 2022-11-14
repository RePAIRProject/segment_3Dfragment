#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <string>
#include <Fragment.h>
#include <Segment.h>



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

	std::string fragmentPath = all_args[0];
	std::string outFileName = all_args[1];

	
	ObjFragment fragment = ObjFragment(fragmentPath);
	std::cout << "Loading data" << std::endl;
	fragment.load();

	std::cout << "Start to segment" << std::endl;
	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;

	bool isSegmented = false;

	double fracture = 0.65;
	std::vector<Segment> segments;
	int nTrials = 1;

	while (!isSegmented)
	{
		double simThresh = fragment.getSimilarThreshByPos(fracture);
		std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
		fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);
		fragment.filterSmallRegions(segments, oRegionsList);
		std::cout << "In trial: " << nTrials <<" Found " << segments.size() << " Segments" << std::endl;

		if (segments.size() == 0)
		{
			fracture = fracture + 0.05;
		}
		else {
			if (segments.size() == 1)
			{
				fracture = fracture - 0.12;
			}
			else
			{
				isSegmented = true;
			}
			
		}

		nTrials += 1;

		if (fracture > 1)
		{
			std::cout << "fracrure value is " << fracture << " and it should be between 0 to 1" << std::endl;
			std::cout << "Please rerun mannually..exiting" << std::endl;
			std::exit(1);
		}

	}

	Segment intactSurface;
	fragment.extractIntactSurface(intactSurface,segments);
	
	std::string fragFolderPath = fragmentPath.substr(0, fragmentPath.find_last_of("\\/"));
	fragment.saveAsObj(fragFolderPath + "\\" + outFileName, intactSurface);
	std::cout << "Write successfully the output to path " << fragFolderPath << std::endl;

	/*
	// For debug:
	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(intactSurface.m_Vertices, intactSurface.m_Faces);
	viewer.data().set_face_based(true);
	viewer.launch();*/
}
