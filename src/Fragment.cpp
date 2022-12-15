#include "Fragment.h"
#include <igl/readOBJ.h>
#include <igl/principal_curvature.h>
#include <igl/adjacency_list.h>
#include <set>
#include <igl/writeOBJ.h>

#include <iostream>
#include <cstdio>
#include <filesystem>

// to remove tmp files

ObjFragment::ObjFragment()
{
}


ObjFragment::ObjFragment(std::string filePath)
{
	m_filePath = filePath;
}


void ObjFragment::load()
{
	Eigen::MatrixXd V;
	igl::readOBJ(m_filePath, V, m_TextureCoordinates, 
		m_Normals, m_Faces, m_Faces2TextureCoordinates, m_Faces2Normals);
	m_FolderPath = m_filePath.substr(0, m_filePath.find_last_of("/\\"));


	m_Vertices = V.block(0, 0, V.rows(), 3);
	m_Colors = V.block(0, 3, V.rows(), 3); // when playing the cub file it is seems to be meaningless
	igl::adjacency_list(m_Faces, m_adjacentVertices); //list of lists containing at index i the adjacent vertices of vertex i
	std::vector < std::vector<int>> _;
	igl::vertex_triangle_adjacency(m_Vertices, m_Faces, m_VerticesAdjacentFacesList, _);


	/*Eigen::MatrixXd pd1, pd2;
	Eigen::VectorXd pv1, pv2;
	igl::principal_curvature(m_Vertices, m_Faces, pd1, pd2, pv1, pv2);
	m_MeshCurvedness = 0.5*(pv1.array().square() + pv2.array().square()).sqrt();*/

	/*std::ofstream f("..\\fragments\\group_39\\processed\\RPf_00320_curvedness.txt");
	for (int i = 0; i < m_MeshCurvedness.rows(); ++i) {
		f << m_MeshCurvedness(i) << '\n';
	}*/

	std::cout << "Warning you read the curvedness HARDCODED" << std::endl;
	m_MeshCurvedness.resize(m_Vertices.rows());
	std::ifstream input("..\\fragments\\group_39\\processed\\RPf_00320_curvedness.txt");
	int i = 0;
	for (std::string line; getline(input, line); )
	{
		m_MeshCurvedness(i) = std::stod(line);
		++i;
	}


	Eigen::MatrixXd curvednessLog = m_MeshCurvedness.array().log();
	//m_NormedMeshCurvedness = 
	double minCurvedness = curvednessLog.minCoeff();// This might be -inf 
	double maxCurvedness = curvednessLog.maxCoeff();
	double inf = std::numeric_limits<double>::infinity();

	if (minCurvedness == -inf || maxCurvedness == inf)
	{
		std::cout << "Min or Max curvdeness are infinity after applying log...cancel log operation" << std::endl;
		//m_NormedMeshCurvedness = m_MeshCurvedness;
		minCurvedness = m_MeshCurvedness.minCoeff();
		maxCurvedness = m_MeshCurvedness.maxCoeff();
		m_NormedMeshCurvedness = ((m_MeshCurvedness.array() - minCurvedness) / (maxCurvedness - minCurvedness));
	}
	else
	{
		m_NormedMeshCurvedness = ((curvednessLog.array() - minCurvedness) / (maxCurvedness - minCurvedness));
	}
	
	


	
}


double ObjFragment::getSimilarThreshByPos(double fracture)
{
	/*
	*  fracture - number between 0 to 1
	*/

	
	std::vector<double> log_norm_curvedness;
	log_norm_curvedness.resize(m_NormedMeshCurvedness.size());
	Eigen::VectorXd::Map(&log_norm_curvedness[0], log_norm_curvedness.size()) = m_NormedMeshCurvedness;
	std::sort(log_norm_curvedness.begin(), log_norm_curvedness.end());
	const auto median_it = log_norm_curvedness.begin() + log_norm_curvedness.size() * fracture;//0.65 good for fine segmentation  //0.5 good for detecting intact vs fractured when removing region merging 
	std::nth_element(log_norm_curvedness.begin(), median_it, log_norm_curvedness.end());
	return *median_it;
}


void ObjFragment::segmentByCurvedness(std::vector<std::vector<int>> &oRegionsList, std::vector<std::vector<int>> &oRegionOutsideBoundaryVerticesList,double similarThreshold)
{

	
	std::map<int, double> available_curves;
	for (auto i = 0; i < m_NormedMeshCurvedness.size(); i++)
	{
		//available_curves[i] = static_cast<float>(mesh_curvedness_[i]);
		available_curves[i] = static_cast<double>(m_NormedMeshCurvedness[i]);
	}

	while (!available_curves.empty())
	{
		std::unordered_map<int, int> current_region;
		std::unordered_map<int, int> current_region_boundary_neighbors;
		std::vector<int> current_seeds;
		current_seeds.reserve(available_curves.size());
		int min_curvature_index;
		//initialize_current_region(available_curves, current_region, current_seeds, min_curvature_index);
		min_curvature_index = available_curves.begin()->first;

		available_curves.erase(min_curvature_index);
		current_region.emplace(min_curvature_index, min_curvature_index); // sort of inserting
		current_seeds.push_back(min_curvature_index);

		grow_current_region(available_curves, current_region, current_region_boundary_neighbors, current_seeds, min_curvature_index, similarThreshold);// mani

		std::vector<int> vec_current_region;
		for (const auto& curr : current_region)
			vec_current_region.push_back(curr.second);
		oRegionsList.push_back(vec_current_region);

		std::vector<int> vec_current_region_boundary_neighbors;
		for (const auto& curr : current_region_boundary_neighbors)
			vec_current_region_boundary_neighbors.push_back(curr.second);
		oRegionOutsideBoundaryVerticesList.push_back(vec_current_region_boundary_neighbors);
	}
}


//int ObjFragment::findIntactSegmentIndex(std::vector<Segment>& segments)
//{
//	int chosenIndex = 0;
//	double minMeanCurvedness = 9999999;
//	int k = 0;
//
//	for (Segment& seg : segments)
//	{
//		double segCurvedness = 0;
//		for (int verIndex : seg.piece_vertices_index_)
//		{
//			segCurvedness += m_NormedMeshCurvedness[verIndex];
//		}
//
//
//		double segAvgCur = segCurvedness / seg.piece_vertices_index_.size();
//
//		if (segAvgCur < minMeanCurvedness)
//		{
//			minMeanCurvedness = segAvgCur;
//			chosenIndex = k;
//		}
//
//		++k;
//	}
//
//	return chosenIndex;
//}


void ObjFragment::grow_current_region(std::map<int, double>& available_curves, std::unordered_map<int, int>& current_region, std::unordered_map<int, int>& current_region_boundary_neighbors, std::vector<int> current_seeds, int min_curvature_index, double segment_threshold_value)
{

	auto segment_avg = m_NormedMeshCurvedness(min_curvature_index);
	auto count = 1;
	for (auto i = 0; i < current_seeds.size(); i++)
	{
		for (auto j = 0; j < m_adjacentVertices[current_seeds[i]].size(); j++)
		{
			auto current_neighbor = m_adjacentVertices[current_seeds[i]][j];
			if ((available_curves.count(current_neighbor) == 1) && (m_NormedMeshCurvedness(current_neighbor) * m_NormedMeshCurvedness(current_seeds[i]) < segment_threshold_value))
			{
				current_region.emplace(current_neighbor, current_neighbor);
				available_curves.erase(current_neighbor);
				if (m_NormedMeshCurvedness(current_neighbor) < segment_threshold_value)
				{
					current_seeds.push_back(current_neighbor);
				}
				else
				{
					//add the neighbors to the region boundary
					for (auto k = 0; k < m_adjacentVertices[current_neighbor].size(); k++)
					{
						auto boundary_neighbor = m_adjacentVertices[current_neighbor][k];
						if (current_region.count(boundary_neighbor) == 0)
						{
							current_region_boundary_neighbors.emplace(boundary_neighbor, boundary_neighbor);
						}
					}

				}
			}
			else
			{
				if (current_region.count(current_neighbor) == 0)
				{
					current_region_boundary_neighbors.emplace(current_neighbor, current_neighbor);
				}
			}

		}
	}
}




//// consider to move this code to the script (segment_intact_surface function)
//void ObjFragment::filterSmallRegions(std::vector<Segment> &segments,std::vector<std::vector<int>> &regions_list_)
//{
//
//	std::map<int, std::vector<int>> small_regions;
//	double minimal_segment_area_percent = 0.05;
//
//	for (auto i = 0; i < regions_list_.size(); i++)
//	{
//		if ((static_cast<double>(regions_list_[i].size()) / static_cast<double>(m_Vertices.rows())) <= minimal_segment_area_percent)
//		{
//			small_regions[i] = regions_list_[i];
//		}
//		else
//		{
//			const Segment current = Segment(regions_list_[i]);
//			segments.push_back(current);
//			//segmented_regions_.push_back(current);
//		}
//	}
//}



//
//void ObjFragment::saveAsObj(std::string outputPath, Segment &iSegment)
//{
//	std::string tmpFilePath = m_filePath + ".tmp"; //"..\\fragments\\cube\\cube_igl_tmp.obj";
//	igl::writeOBJ(tmpFilePath, iSegment.m_Vertices, iSegment.m_Faces, m_Normals,
//		iSegment.m_Faces2Normals, m_TextureCoordinates, iSegment.m_Faces2TextureCoordinates);
//
//	std::ifstream input(m_filePath);
//	std::string mtlibLine;
//	std::string materialLine;
//
//	for (std::string line; getline(input, line); )
//	{
//		if (line.find("mtllib") != std::string::npos)
//		{
//			mtlibLine = line;
//		}
//
//		if (line.find("usemtl") != std::string::npos)
//		{
//			materialLine = line;
//		}
//
//		if (mtlibLine.length() > 0 && materialLine.length() > 0)
//		{
//			break;
//		}
//
//	}
//
//	std::ofstream outFile(outputPath);
//	outFile << "# Generated with Libigl \n" << mtlibLine << '\n' << materialLine << "\n\n";
//
//	{
//		std::ifstream objData(tmpFilePath);
//
//		for (std::string line; getline(objData, line); )
//		{
//			outFile << line << "\n";
//		}
//	}
//
//	std::filesystem::remove(tmpFilePath);
//}


double ObjFragment::localAvgCurvedness(const std::vector<int>& vertIndexes)
{
	double sum = 0;
	for (int i = 0; i < vertIndexes.size(); i++)
	{
		sum += m_MeshCurvedness[vertIndexes[i]];
	}

	return sum / vertIndexes.size();
}

Eigen::Vector3d ObjFragment::localAvgNormal(const std::vector<int>& vertIndexes)
{
	Eigen::Vector3d sum = Eigen::Vector3d::Zero();
	for (auto iVert: vertIndexes)
	{
		for (int i = 0; i <3; i++)
		{
			sum[i] += m_Normals.coeff(iVert, i);
		}
	}

	return sum / vertIndexes.size();
}
