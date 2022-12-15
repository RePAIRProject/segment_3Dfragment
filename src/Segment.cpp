#include <Segment.h>

Segment::Segment()
{

}

Segment::Segment(std::vector<int> piece_vertices_index, std::vector<int> piece_vertices_indexes_boundary)
{
	piece_vertices_index_ = piece_vertices_index;
	m_OutsideBoundaryVertsIndexes = piece_vertices_indexes_boundary;
}


void Segment::loadBasicData( const std::vector <std::vector<int>> &parentVerticesAdjacentFacesList,
	 const Eigen::MatrixXi &parentFaces, const Eigen::MatrixXd &parentVertices,
	 const Eigen::MatrixXi &parentFaces2TextureCoordinates,
	 const Eigen::MatrixXi &parentFaces2Normals)
{

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
	std::vector<int>FragVertIndexes_;
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


