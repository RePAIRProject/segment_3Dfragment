#include <Segment.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <cstdio>
#include <filesystem>

ObjFragment emptyFrag;
Segment::Segment():m_ParentFragment(emptyFrag)
{

}

Segment::Segment(std::vector<int> piece_vertices_index, 
			std::vector<int> piece_vertices_indexes_boundary,  ObjFragment& parentFragment) : m_ParentFragment(parentFragment)
{
	piece_vertices_index_ = piece_vertices_index;
	m_OutsideBoundaryVertsIndexes = piece_vertices_indexes_boundary;
}


double Segment::calcAvgCurvedness()
{
	m_AvgCurvedness =  ObjFragment::localAvgCurvedness(piece_vertices_index_);
	return m_AvgCurvedness;
}


/// This need to be fixed
void Segment::loadNormedNormals()
{
	//m_NormedNormals = Eigen::MatrixXd::Zero(piece_vertices_index_.size(), 3);
	const Eigen::MatrixXd& parentNormals = m_ParentFragment.m_Normals;
	for (int vertexIndex : piece_vertices_index_)
	{
		Eigen::Vector3d norm = { parentNormals(vertexIndex, 0),
									parentNormals(vertexIndex, 1),
										parentNormals(vertexIndex, 2) };

		m_NormedNormals.push_back(norm);

	}
}


void Segment::loadBasicData()
{
	const std::vector <std::vector<int>>& parentVerticesAdjacentFacesList = m_ParentFragment.m_VerticesAdjacentFacesList;
	const Eigen::MatrixXi& parentFaces = m_ParentFragment.m_Faces;
	const Eigen::MatrixXd& parentVertices = m_ParentFragment.m_Vertices;
	const Eigen::MatrixXi& parentFaces2TextureCoordinates = m_ParentFragment.m_Faces2TextureCoordinates;
	const Eigen::MatrixXi& parentFaces2Normals = m_ParentFragment.m_Faces2Normals;
	const Eigen::MatrixXd& parentNormals = m_ParentFragment.m_Normals;
	m_ParentFilePath = m_ParentFragment.m_filePath;

	std::set<int> uniqueVerticesIndexes;
	std::set<int> uniqueFacesIndexes;
	for (int vertexIndex : piece_vertices_index_)
	{
		//int  = oIntactSurface.piece_vertices_index_[i];
		//uniqueVerticesIndexes.insert(vertexIndex);
		for (int faceIndex : parentVerticesAdjacentFacesList[vertexIndex])
		{
			uniqueFacesIndexes.insert(faceIndex);

			for (auto j = 0; j < parentFaces.cols(); j++)
			{
				uniqueVerticesIndexes.insert(parentFaces.coeff(faceIndex, j));
			}
		}
	}

	std::vector<int> FragfacesIndexes_;
	std::copy(uniqueFacesIndexes.begin(), uniqueFacesIndexes.end(), std::back_inserter(FragfacesIndexes_));
	std::vector<int>FragVertIndexes_; // note: this is not equal size to piece_vertices_index_ (even larger!)
	std::copy(uniqueVerticesIndexes.begin(), uniqueVerticesIndexes.end(), std::back_inserter(FragVertIndexes_));

	int nSize = FragVertIndexes_.size();
	m_Vertices.resize(nSize, 3);
	std::map<int, int> VertIndex2SegIndex;

	int i = 0;
	for (int vertexIndex : FragVertIndexes_)
	{
		VertIndex2SegIndex.insert({ vertexIndex,i });
		for (int j = 0; j < 3; ++j)
			m_Vertices(i, j) = parentVertices(vertexIndex, j);

		i++;
	}

	m_Normals = Eigen::MatrixXd::Zero(nSize, 3); // .resize(nSize, 3);

	i = 0;
	for (int vertexIndex : FragVertIndexes_)
	{
		for (int j = 0; j < 3; ++j)
			m_Normals(i, j) = parentNormals(vertexIndex, j);

		i++;
	}


	m_Faces.resize(FragfacesIndexes_.size(), 3);
	//Eigen::MatrixXi FTC, FN;
	m_Faces2Normals.resize(FragfacesIndexes_.size(), 3);
	m_Faces2TextureCoordinates.resize(FragfacesIndexes_.size(), 3);

	i = 0;
	for (int ixFragfaces : FragfacesIndexes_)
	{
		for (int j = 0; j < 3; ++j)
		{
			int fragVertIndex = parentFaces(ixFragfaces, j);
			m_Faces(i, j) = VertIndex2SegIndex[fragVertIndex];

			fragVertIndex = parentFaces2TextureCoordinates(ixFragfaces, j);
			m_Faces2TextureCoordinates(i, j) = fragVertIndex; //VertIndex2SegIndex[fragVertIndex];

			fragVertIndex = parentFaces2Normals(ixFragfaces, j);
			m_Faces2Normals(i, j) = fragVertIndex; //VertIndex2SegIndex[fragVertIndex];

		}
		++i;
	}
}

void Segment::saveAsObj(std::string outputPath)
{
	auto &parentNormals = m_ParentFragment.m_Normals;
	const Eigen::MatrixXd& parentTextureCoordinates = m_ParentFragment.m_TextureCoordinates;

	std::string tmpFilePath = m_ParentFilePath + ".tmp"; //"..\\fragments\\cube\\cube_igl_tmp.obj";
	igl::writeOBJ(tmpFilePath, m_Vertices, m_Faces, m_Normals, //parentNormals,
		m_Faces2Normals, parentTextureCoordinates, m_Faces2TextureCoordinates);

	std::ifstream input(m_ParentFilePath);
	std::string mtlibLine;
	std::string materialLine;

	for (std::string line; getline(input, line); )
	{
		if (line.find("mtllib") != std::string::npos)
		{
			mtlibLine = line;
		}

		if (line.find("usemtl") != std::string::npos)
		{
			materialLine = line;
		}

		if (mtlibLine.length() > 0 && materialLine.length() > 0)
		{
			break;
		}

	}

	std::ofstream outFile(outputPath);
	outFile << "# Generated with Libigl \n" << mtlibLine << '\n' << materialLine << "\n\n";

	{
		std::ifstream objData(tmpFilePath);

		for (std::string line; getline(objData, line); )
		{
			outFile << line << "\n";
		}
	}

	std::filesystem::remove(tmpFilePath);
}
