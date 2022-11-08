#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <string>
#include <Fragment.h>
#include <Segment.h>



int main(int argc, char* argv[])
{
	std::string fragmentName = "RPf_00154"; //"cube";

	std::string fragmentPath = "..\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";//".\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";
	std::string outputPath = "..\\fragments\\" + fragmentName + "\\" ;
	ObjFragment fragment = ObjFragment(fragmentPath);
	fragment.load();

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;

	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, 0.1);
	std::vector<Segment> segments;
	fragment.filterSmallRegions(segments, oRegionsList);

	Segment intactSurface;
	Eigen::VectorXd normedCurvedness;//.log();
	normedCurvedness = fragment.m_MeshCurvedness.array().log();
	normedCurvedness = ((normedCurvedness.array() - normedCurvedness.minCoeff()) / (normedCurvedness.maxCoeff() - normedCurvedness.minCoeff()));
	fragment.extractIntactSurface(intactSurface,segments, normedCurvedness);
	
	fragment.saveAsObj(outputPath + fragmentName + "_igl.obj", intactSurface);

	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(intactSurface.m_Vertices, intactSurface.m_Faces);
	viewer.data().set_face_based(true);
	viewer.launch();
}
