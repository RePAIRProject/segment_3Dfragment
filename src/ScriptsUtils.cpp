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



void mergeSmall2BigSegments(std::map<int, Segment*> &smallSegments,std::map<int, Segment*> &bigSegments, 
							std::map<int, int>& vertIndex2SegIndex)
{
	int nLastFreeDebug = 0;
	std::set<int> ixCurrSegBigNeigh;
	std::set<int> ixCurrSegSmallNeigh;
	std::map<int, int> segSrc2segDst;

	while (!smallSegments.empty()) //debugIter++ < 1()
	{
		
		if (nLastFreeDebug == smallSegments.size())
		{
			std::cout << "WARNING: mergeSmall2BigSegments function did not converged after " << std::endl;
			break;
		}
		nLastFreeDebug = smallSegments.size();

		//std::cout << "WARNING: merge only for a single big segment" << std::endl;
		//for (auto bigSegIt = std::next(bigSegments.begin(), ixMaxSize_debug); bigSegIt != std::next(bigSegments.begin(), ixMaxSize_debug+1); bigSegIt++)//bigSegments.end()
		for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)//bigSegments.end()
		{			
			/*
				Compute the neighboors segments
			*/
			for (int iVertBoundary : bigSegIt->second->m_OutsideBoundaryVertsIndexes)
			{
				int iNeighboorSeg = vertIndex2SegIndex[iVertBoundary];

				if (bigSegments.count(iNeighboorSeg) > 0)
				{
					ixCurrSegBigNeigh.insert(iNeighboorSeg);
				}
				else
				{
					ixCurrSegSmallNeigh.insert(iNeighboorSeg);
				}
			}

			for (std::set<int>::iterator itSmallNeigh = ixCurrSegSmallNeigh.begin(); itSmallNeigh != ixCurrSegSmallNeigh.end(); itSmallNeigh++)
			{
				segSrc2segDst[*itSmallNeigh] = bigSegIt->first;
			}
			

			ixCurrSegBigNeigh.clear();
			ixCurrSegSmallNeigh.clear();
		}

		/*
			Merging
		*/
		for (std::map<int,int>::iterator itSegSrc2Dst = segSrc2segDst.begin(); itSegSrc2Dst !=segSrc2segDst.end();++itSegSrc2Dst)
		{
			int iSrcSeg = itSegSrc2Dst->first;
			int iDstSeg = itSegSrc2Dst->second;

			
			if (smallSegments.count(iSrcSeg) > 0)
			{
				merge(iSrcSeg, iDstSeg, vertIndex2SegIndex, smallSegments, bigSegments);
			}
		}

		for (std::map<int, int>::iterator itSrc2Dst = segSrc2segDst.begin(); itSrc2Dst != segSrc2segDst.end(); ++itSrc2Dst)
		{
			smallSegments.erase(itSrc2Dst->first);
		}
		segSrc2segDst.clear();
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
