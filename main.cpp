#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
#include <string>
#include <Fragment.h>
#include <Segment.h>



int main(int argc, char* argv[])
{
	std::string fragmentName = "RPf_00154"; //"cube";

	//std::string fragmentPath = "..\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";//".\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";
	std::string fragmentPath = "..\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";//".\\fragments\\" + fragmentName + "\\" + fragmentName + ".obj";
	std::string outputPath = "..\\fragments\\" + fragmentName + "\\" ;
	ObjFragment fragment = ObjFragment(fragmentPath);
	fragment.load();
	std::vector<Eigen::MatrixXd> Vs;
	std::vector<Eigen::MatrixXd> Fs;

	Segment s;
	fragment.extractIntactSurface(Vs, Fs,s);
	fragment.saveAsObj(outputPath + fragmentName + "_igl.obj",s);

	igl::opengl::glfw::Viewer viewer;
	viewer.data().set_mesh(s.m_Vertices, s.m_Faces);
	viewer.data().set_face_based(true);
	viewer.launch();
}
