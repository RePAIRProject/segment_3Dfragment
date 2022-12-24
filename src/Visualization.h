#pragma once
#include <igl/opengl/glfw/Viewer.h>


class Visualizer
{
public:
	Visualizer();
	int appendMesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& colors);
	void launch();
	void writeOFF(std::string fileName, Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd C);
//private:
	igl::opengl::glfw::Viewer m_Viewer;
};
void generateRandomColors(std::map<int, Eigen::RowVector3d> &oColors, int numberColors);