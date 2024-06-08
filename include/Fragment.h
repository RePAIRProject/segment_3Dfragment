#pragma once
#include <string>
#include <igl/principal_curvature.h>
//#include <../src/Segment.h>
// Todo : make an fragment class that will have V and F only


class ObjFragment
{

public:
	ObjFragment();
	ObjFragment(std::string filePath);
	
	std::string m_filePath;
	std::string m_FolderPath;
	std::string m_OutputPath;
	std::string m_Name;

	Eigen::MatrixXd m_Vertices;//   V  double matrix of vertex positions  #V by 3
	Eigen::MatrixXd m_TextureCoordinates;//   TC  double matrix of texture coordinats #TC by 2
	Eigen::MatrixXd m_Normals;//   N  double matrix of corner normals #N by 3
	Eigen::MatrixXi m_Faces;//   F  #F list of face indices into vertex positions
	Eigen::MatrixXi m_Faces2TextureCoordinates;//   FTC  #F list of face indices into vertex texture coordinates
	Eigen::MatrixXi m_Faces2Normals;//   FN  #F list of face indices into vertex normals
	Eigen::MatrixXd m_Colors;
	Eigen::VectorXd m_MeshCurvedness;
	Eigen::VectorXd m_NormedMeshCurvedness;
	std::vector<std::vector<int>> m_adjacentVertices;
	std::vector <std::vector<int>> m_VerticesAdjacentFacesList; // indexes from zero...

	void load();
	
	double getSimilarThreshByPos(double fracture);
	void segmentByCurvedness(std::vector<std::vector<int>>& oRegionsList,std::vector<std::vector<int>>& oRegionOutsideBoundaryVerticesList, double similarThreshold);
	void grow_current_region(std::map<int, double>& available_curves, std::unordered_map<int, int>& current_region, std::unordered_map<int, int>& current_region_boundary_neighbors, std::vector<int> current_seeds, int min_curvature_index, double segment_threshold_value);

	double localAvgCurvedness(const std::vector<int> &vertIndexes);
	Eigen::Vector3d localAvgNormal(const std::vector<int>& vertIndexes);
	
	void crop(std::vector<int> vertRemainIndexes);


};



