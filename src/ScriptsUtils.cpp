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

void sortToSmallAndBigSegments(std::map<int, Segment*>& oSmallSegments, std::map<int, Segment*>& oBigSegments,
								std::vector<Segment> &segments, double minBigSegPercSize)
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

void merge(int iSrcSeg, int iDstSeg, std::map<int, int> &vertIndex2SegIndex,
			std::map<int, Segment*> &srcSegPool, std::map<int, Segment*> &dstSegPool)
{
	for (int iVert : srcSegPool.at(iSrcSeg)->piece_vertices_index_)
	{
		vertIndex2SegIndex[iVert] = iDstSeg;
	}

	dstSegPool.at(iDstSeg)->piece_vertices_index_.insert(
		dstSegPool.at(iDstSeg)->piece_vertices_index_.end(),
		srcSegPool.at(iSrcSeg)->piece_vertices_index_.begin(),
		srcSegPool.at(iSrcSeg)->piece_vertices_index_.end()
	);
		
	dstSegPool.at(iDstSeg)->m_OutsideBoundaryVertsIndexes.insert(
		dstSegPool.at(iDstSeg)->m_OutsideBoundaryVertsIndexes.end(),
		srcSegPool.at(iSrcSeg)->m_OutsideBoundaryVertsIndexes.begin(),
		srcSegPool.at(iSrcSeg)->m_OutsideBoundaryVertsIndexes.end()
	);
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

Eigen::Vector3d calcVariance(const std::vector<Eigen::Vector3d>& vectors, Eigen::Vector3d mean)
{
	Eigen::Vector3d sum = Eigen::Vector3d::Zero();
	for (auto& vecIt : vectors)
	{
		for (int i = 0; i < 3; i++)
		{
			sum[i] = sum[i] + (vecIt[i] - mean[i])* (vecIt[i] - mean[i]);
		}
	}

	return sum / vectors.size();
}
