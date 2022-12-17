#pragma once
#include <Segment.h>
#include <Fragment.h>

void initSegments(std::vector<Segment> &oSegments, std::map<int, int> &oVertIndex2SegIndex,
					const std::vector<std::vector<int>> &oRegionsList,
					const std::vector<std::vector<int>> &oRegionOutsideBoundaryVerticesList,
					 double fragmentSize, ObjFragment& parentFragment);

void sortToSmallAndBigSegments(std::map<int, Segment*>& oSmallSegments, std::map<int, Segment*>& oBigSegments, 
							std::vector<Segment> &segments, double minBigSegPercSize);

Eigen::Vector3d calcAvg(const std::vector<Eigen::Vector3d>& vectors);
Eigen::Vector3d calcVariance(const std::vector<Eigen::Vector3d>& vectors, Eigen::Vector3d mean);
void merge(int iSrcSeg, int iDstSeg, std::map<int, int>& vertIndex2SegIndex, std::map<int, Segment*>& srcSegPool, std::map<int, Segment*>& dstSegPool);
void mergeSmall2BigSegments(std::map<int, Segment*>& smallSegments, std::map<int, Segment*>& bigSegments, std::map<int, int>& vertIndex2SegIndex);
