#pragma once
#include <igl/opengl/glfw/Viewer.h>


class Visualizer
{
public:
	Visualizer();
	void generateRandomColors(std::map<int, Eigen::RowVector3d> &oColors, int numberColors);
	int appendMesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& colors);
	void launch();
private:
	igl::opengl::glfw::Viewer m_Viewer;
};