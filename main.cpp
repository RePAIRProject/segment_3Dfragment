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

	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, 0.1);
	std::vector<Segment> segments;
	fragment.filterSmallRegions(segments, oRegionsList);

	Segment intactSurface;
	Eigen::VectorXd normedCurvedness;
	normedCurvedness = fragment.m_MeshCurvedness.array().log();
	normedCurvedness = ((normedCurvedness.array() - normedCurvedness.minCoeff()) / (normedCurvedness.maxCoeff() - normedCurvedness.minCoeff()));
	fragment.extractIntactSurface(intactSurface,segments, normedCurvedness);
	
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
