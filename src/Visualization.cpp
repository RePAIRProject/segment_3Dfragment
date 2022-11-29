#include <Visualization.h>

Visualizer::Visualizer()
{
	
}

void Visualizer::generateRandomColors(std::map<int, Eigen::RowVector3d> &oColors,int numberColors)
{
	for (int i = 0; i < numberColors; i++)
	{
		oColors.emplace(i, 0.5 * Eigen::RowVector3d::Random().array() + 0.5);
	}
}

int Visualizer::appendMesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& colors)
{

	//std::map<int, Eigen::RowVector3d> colors;
	//int segID = viewer.append_mesh();
	//viewer.data(segID).set_mesh(intactSurface.m_Vertices, intactSurface.m_Faces);
	//colors.emplace(viewer.data().id, 0.5 * Eigen::RowVector3d::Random().array() + 0.5);
	//viewer.data(segID).set_colors(colors[segID]);
	//viewer.data(segID).show_lines = false;

	//viewer.data().set_face_based(true);
	//viewer.launch();

	
	int segID = m_Viewer.append_mesh();
	m_Viewer.data(segID).set_mesh(vertices, faces);
	m_Viewer.data(segID).set_colors(colors);
	m_Viewer.data(segID).show_lines = false;

	return segID;
}

void Visualizer::launch()
{
	m_Viewer.data().set_face_based(true);
	m_Viewer.launch();
}