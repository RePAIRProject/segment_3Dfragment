#pragma once
#include <Segment.h>


void initSegments(std::vector<Segment> &oSegments, std::map<int, int> &oVertIndex2SegIndex,
					const std::vector<std::vector<int>> &oRegionsList,
					const std::vector<std::vector<int>> &oRegionOutsideBoundaryVerticesList,
					 double fragmentSize);

void sortToSmallAndBigSegments(std::map<int, Segment*>& oSmallSegments, std::map<int, Segment*>& oBigSegments, std::vector<Segment> segments, double minBigSegPercSize);

