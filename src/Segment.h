#pragma once
#include <igl/adjacency_list.h>


class Segment
{
public:
	Segment();

	std::vector<int> piece_vertices_index_;
	
	Segment(std::vector<int> vertices_index);

	Segment(Eigen::MatrixXd V, Eigen::MatrixXd  TC, Eigen::MatrixXd N,
		Eigen::MatrixXi F, Eigen::MatrixXi  FTC, Eigen::MatrixXi FN);

	Eigen::MatrixXd m_Vertices;//   V  double matrix of vertex positions  #V by 3
	Eigen::MatrixXd m_TextureCoordinates;//   TC  double matrix of texture coordinats #TC by 2
	Eigen::MatrixXd m_Normals;//   N  double matrix of corner normals #N by 3
	Eigen::MatrixXi m_Faces;//   F  #F list of face indices into vertex positions
	Eigen::MatrixXi m_Faces2TextureCoordinates;//   FTC  #F list of face indices into vertex texture coordinates
	Eigen::MatrixXi m_Faces2Normals;//   FN  #F list of face indices into vertex normals
	Eigen::MatrixXd m_Colors;

};