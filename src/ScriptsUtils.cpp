#include <ScriptsUtils.h>
#include <Segment.h>

void initSegments(std::vector<Segment>& oSegments, std::map<int, int>& oVertIndex2SegIndex,
	const std::vector<std::vector<int>>& oRegionsList,
	const std::vector<std::vector<int>>& oRegionOutsideBoundaryVerticesList,
	double fragmentSize, ObjFragment& parentFragment)
{
	for (int iSeg = 0; iSeg < oRegionsList.size(); iSeg++)
	{
		Segment currSeg(oRegionsList[iSeg], oRegionOutsideBoundaryVerticesList[iSeg],parentFragment);
		currSeg.m_fracSizeOfFragment = static_cast<double>(currSeg.piece_vertices_index_.size())/fragmentSize;

		for (int jVert = 0; jVert < oRegionsList[iSeg].size(); jVert++)
		{
			oVertIndex2SegIndex.insert({ oRegionsList[iSeg][jVert],iSeg });
		}

		oSegments.push_back(currSeg);
	}
}

void sortToSmallAndBigSegments(std::map<int, Segment*>& oSmallSegments, std::map<int, Segment*>& oBigSegments, std::vector<Segment> segments, double minBigSegPercSize)
{
	for (int iSeg = 0; iSeg < segments.size(); iSeg++)
	{
		
		if (segments[iSeg].m_fracSizeOfFragment < minBigSegPercSize) {
			oSmallSegments.insert({ iSeg, &segments[iSeg] });
		}
		else
		{
			oBigSegments.insert({ iSeg, &segments[iSeg] });
		}
		
	}
}

Eigen::Vector3d calcAvg(const std::vector<Eigen::Vector3d>& vectors)
{

	Eigen::Vector3d sum = Eigen::Vector3d::Zero();
	for (auto &vecIt : vectors)
	{
		sum = sum + vecIt;
	}

	return sum / vectors.size();
}
