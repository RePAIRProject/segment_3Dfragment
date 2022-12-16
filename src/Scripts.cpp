#include "Scripts.h"
#include "Visualization.h"
#include "ScriptsUtils.h"

void segment_intact_surface(std::vector<std::string> all_args)
{
	std::string fragmentPath = all_args[0];
	std::string outFileName = all_args[1];


	ObjFragment fragment = ObjFragment(fragmentPath);
	std::cout << "Loading data" << std::endl;
	fragment.load();

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	std::map<int, Segment*> smallSegments;
	std::map<int, Segment*> bigSegments;

	bool isSegmented = false;
	int intactIndex = -1;
	double fracture = 0.65;
	double minBigSegPercSize = 0.05;
	int nTrials = 1;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	while (!isSegmented)
	{
		oRegionsList.clear();
		oRegionOutsideBoundaryVerticesList.clear();
		vertIndex2SegIndex.clear();
		segments.clear();
		smallSegments.clear();
		bigSegments.clear();

		double simThresh = fragment.getSimilarThreshByPos(fracture);
		std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
		fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);

		initSegments(segments, vertIndex2SegIndex, oRegionsList, oRegionOutsideBoundaryVerticesList, fragmentSize, fragment);
		sortToSmallAndBigSegments(smallSegments, bigSegments, segments, minBigSegPercSize);

		//fragment.filterSmallRegions(segments, oRegionsList);
		std::cout << "In trial: " << nTrials << " Found " << bigSegments.size() << "big Segments" << std::endl;

		if (bigSegments.size() == 0)
		{
			fracture = fracture + 0.05;
			continue;
		}
		else {
			if (bigSegments.size() == 1)
			{
				fracture = fracture - 0.12;
				continue;
			}
			else
			{
				//isSegmented = true;
			}

		}

		nTrials += 1;

		if (fracture > 1)
		{
			std::cout << "fracrure value is " << fracture << " and it should be between 0 to 1" << std::endl;
			std::cout << "Please rerun mannually..exiting" << std::endl;
			std::exit(1);
		}

		intactIndex = -1;
		double minMeanCurvedness = 9999999;
		int k = 0;

		for (Segment& seg : segments)
		{
			double segCurvedness = 0;
			for (int verIndex : seg.piece_vertices_index_)
			{
				segCurvedness += fragment.m_NormedMeshCurvedness[verIndex];
			}

			double segAvgCur = segCurvedness / seg.piece_vertices_index_.size();

			if (segAvgCur < minMeanCurvedness)
			{
				minMeanCurvedness = segAvgCur;
				intactIndex = k;
			}

			++k;
		}

		segments[intactIndex].loadNormedNormals();
		Eigen::Vector3d avgNormal = calcAvg(segments[intactIndex].m_NormedNormals);
		Eigen::Vector3d stdNormal = calcVariance(segments[intactIndex].m_NormedNormals, avgNormal).array().sqrt();
		double l2 = stdNormal.norm();

		if (l2 < 0.2)
		{
			isSegmented = true;
		}
		else {
			fracture = fracture - 0.07;
		}
	}

	

	segments[intactIndex].loadBasicData();
	segments[intactIndex].saveAsObj(fragment.m_FolderPath + "\\" + outFileName);
	std::cout << "Write successfully the output to path " << fragment.m_FolderPath << std::endl;

}

//void segment(std::vector<std::string> all_args)
//{
//	std::string fragmentPath = all_args[0];
//	std::string outFileName = all_args[1];
//	ObjFragment fragment = ObjFragment(fragmentPath);
//	std::cout << "Loading data" << std::endl;
//	fragment.load();
//
//	std::cout << "Start to segment" << std::endl;
//	double fracture = 0.45; //0.45; //0.5
//	double simThresh = fragment.getSimilarThreshByPos(fracture);
//	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
//	std::vector<std::vector<int>> oRegionsList;
//	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
//	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);
//
//
//	/*
//	
//		Start mergnig algo - you can put it in other function when refactor
//	*/
//
//	/*
//	* init data structures
//	*/
//	double minSegPercSize = 0.00025;  //0.000125;//0.005 // This need to get as a parameter to the script
//	std::map<int, Segment> smallSegments;
//	std::map<int, Segment> bigSegments;
//	std::map<int, double> segmentsAvgCurvedness; // We avg? maybe std is a parameter we should consider
//	
//	std::map<int, int> vertIndex2SegIndex;
//	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());
//
//	for (int iSeg = 0; iSeg < oRegionsList.size(); iSeg++)
//	{
//		auto segVertsIndexes = oRegionsList[iSeg];
//		Segment currSeg(segVertsIndexes);
//		segmentsAvgCurvedness.insert({ iSeg,fragment.localAvgCurvedness(segVertsIndexes) });
//		
//
//		double segmentSize = static_cast<double>(currSeg.piece_vertices_index_.size());
//		if (segmentSize/fragmentSize < minSegPercSize){
//			smallSegments.insert({ iSeg, currSeg });
//			
//		}
//		else
//		{
//			bigSegments.insert({ iSeg, currSeg });
//
//		}
//
//		for (int i = 0; i < segVertsIndexes.size(); i++)
//		{
//			vertIndex2SegIndex.insert({ segVertsIndexes[i],iSeg });
//		}
//
//	}
//
//	/*
//		For debug
//	*/
//	int ixMaxSize_debug = -1;
//	double maxSizeSegment_debug = 0;
//	int k = 0;
//	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
//	{
//		double segmentSize = static_cast<double>(bigSegIt->second.piece_vertices_index_.size());
//		if (segmentSize > maxSizeSegment_debug)
//		{
//			maxSizeSegment_debug = segmentSize;
//			ixMaxSize_debug = k;
//		}
//		++k;
//	}
//
//	/*
//		Coloring
//	*/
//	Visualizer visualizer;
//	std::map<int, Eigen::RowVector3d> meshColors;
//	visualizer.generateRandomColors(meshColors, oRegionsList.size());
//	visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, meshColors[0]);
//
//	Eigen::MatrixXd bigSegmentColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
//
//	auto colorIt = meshColors.begin();
//	//for (auto bigSegIt = std::next(bigSegments.begin(), ixMaxSize_debug); bigSegIt != std::next(bigSegments.begin(), ixMaxSize_debug + 1); bigSegIt++)
//	for (auto bigSegIt = bigSegments.begin(); bigSegIt !=bigSegments.end(); bigSegIt++)
//	{
//		std::vector<int> regionVerticesIndx = oRegionsList[bigSegIt->first];
//		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//		{
//
//			bigSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//		}
//		++colorIt;
//	}
//
//	Eigen::MatrixXd rawSegmentColors;
//	rawSegmentColors.resize(fragment.m_Vertices.rows(), 4);
//	//segmentBigSmallColors.resize(fragment.m_Vertices.rows(), 4);
//	colorIt = meshColors.begin();
//	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
//	{
//		std::vector<int> regionVerticesIndx = oRegionsList[bigSegIt->first];
//		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//		{
//			rawSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//
//		}
//
//		colorIt++;
//	}
//
//	for (auto smallSegIt = smallSegments.begin(); smallSegIt != smallSegments.end(); smallSegIt++)
//	{
//		std::vector<int> regionVerticesIndx = oRegionsList[smallSegIt->first];
//		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//		{
//			rawSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//		}
//		colorIt++;
//	}
//
//
//	std::set<int> ixCurrSegBigNeigh;
//	std::set<int> ixCurrSegSmallNeigh;
//	std::map<int, int> segSrc2segDst; // The src seg that will be merge to dst seg
//	std::map<int, std::vector<int>> segSrc2MultiDst;
//
//	std::cout << "Start the merge process" << std::endl;
//	int nLastFreeDebug = -1;
//	int debugIter = 0;
//	
//	while (!smallSegments.empty()) //debugIter++ < 1()
//	{
//		debugIter++;
//		if (nLastFreeDebug == smallSegments.size())
//		{
//			std::cout << "WARNING: merge did not converged after " << debugIter << " iterations" << std::endl;
//			break;
//		}
//		nLastFreeDebug = smallSegments.size();
//
//		//std::cout << "WARNING: merge only for a single big segment" << std::endl;
//		//for (auto bigSegIt = std::next(bigSegments.begin(), ixMaxSize_debug); bigSegIt != std::next(bigSegments.begin(), ixMaxSize_debug+1); bigSegIt++)//bigSegments.end()
//		for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)//bigSegments.end()
//		{
//			auto& iVertsAtBoundary = oRegionOutsideBoundaryVerticesList[bigSegIt->first];
//			
//			/*
//				Compute the neighboors segments
//			*/
//			for (int iVertBoundary : iVertsAtBoundary)
//			{
//				int iNeighboorSeg = vertIndex2SegIndex[iVertBoundary];
//
//				// Identify wheter its a big or small segment
//				if (bigSegments.count(iNeighboorSeg) > 0)
//				{
//					ixCurrSegBigNeigh.insert(iNeighboorSeg);
//				}
//				else
//				{
//					ixCurrSegSmallNeigh.insert(iNeighboorSeg);
//				}
//			}
//
//			for (std::set<int>::iterator itSmallNeigh = ixCurrSegSmallNeigh.begin(); itSmallNeigh != ixCurrSegSmallNeigh.end(); itSmallNeigh++)
//			{
//
//				segSrc2segDst[*itSmallNeigh] = bigSegIt->first;
//				//segSrc2segDst.insert({ *itSmallNeigh,bigSegIt->first});
//				//segSrc2MultiDst[*itSmallNeigh].push_back
//			}
//			
//
//			ixCurrSegBigNeigh.clear();
//			ixCurrSegSmallNeigh.clear();
//		}
//
//		//std::cout << "before applying the merging" << std::endl;
//		/*
//			Merging
//		*/
//		for (std::map<int,int>::iterator itSegSrc2Dst = segSrc2segDst.begin(); itSegSrc2Dst !=segSrc2segDst.end();++itSegSrc2Dst)
//		{
//			int iSrcSeg = itSegSrc2Dst->first;
//			int iDstSeg = itSegSrc2Dst->second;
//
//			
//			if (smallSegments.count(iSrcSeg) > 0)
//			{
//
//
//				for (int i = 0; i < smallSegments.at(iSrcSeg).piece_vertices_index_.size(); i++)
//				{
//					vertIndex2SegIndex[smallSegments.at(iSrcSeg).piece_vertices_index_[i]] = iDstSeg;
//				}
//
//				if (bigSegments.count(iDstSeg) > 0)
//				{
//					bigSegments.at(iDstSeg).piece_vertices_index_.insert(
//						bigSegments.at(iDstSeg).piece_vertices_index_.end(),
//						smallSegments.at(iSrcSeg).piece_vertices_index_.begin(),
//						smallSegments.at(iSrcSeg).piece_vertices_index_.end()
//					);
//
//					oRegionOutsideBoundaryVerticesList[iDstSeg].insert(
//						oRegionOutsideBoundaryVerticesList[iDstSeg].end(),
//						oRegionOutsideBoundaryVerticesList[iSrcSeg].begin(),
//						oRegionOutsideBoundaryVerticesList[iSrcSeg].end()
//					);
//				}
//				else {
//					std::cout << "You should not be here";
//				}
//
//			}
//			
//			
//		}
//
//		for (std::map<int, int>::iterator itSrc2Dst = segSrc2segDst.begin(); itSrc2Dst != segSrc2segDst.end(); ++itSrc2Dst)
//		{
//			smallSegments.erase(itSrc2Dst->first);
//		}
//		segSrc2segDst.clear();
//	}
//	
//	/*
//		Coloring
//	*/
//
//	Eigen::MatrixXd segmentColorsAfterMerge = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
//	colorIt = meshColors.begin();
//
//	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
//	{
//		std::vector<int> regionVerticesIndx = bigSegIt->second.piece_vertices_index_;
//		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//		{
//			segmentColorsAfterMerge.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//		}
//		++colorIt;
//	}
//
//	/*
//		Coloring
//	*/
//	Eigen::MatrixXd segmentBoundaryAfterMerge = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
//	colorIt = meshColors.begin();
//
//	for (int iRegion = 0; iRegion < oRegionOutsideBoundaryVerticesList.size(); iRegion++)
//	{
//		if (bigSegments.count(iRegion) > 0)
//		{
//
//			std::vector<int> regionVerticesIndx = oRegionOutsideBoundaryVerticesList[iRegion];
//			for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//			{
//				segmentBoundaryAfterMerge.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//			}
//			++colorIt;
//		}
//	}
//
//	/*
//		Merge the big segments by normals 
//	*/
//
//	if (bigSegments.size() < 6)
//	{
//		std::cout << "ERROR: Expect that the final number of big segment to be greater than 6 (Prior Knoledge)" << std::endl;
//		exit(1);
//	}
//
//	std::map<int, Eigen::Vector3d> segmentsAvgNormal;
//	std::map<int, Eigen::Vector3d> segmentsCenterMass;
//	std::map<int, double> segmentsSize;
//	segmentsAvgCurvedness.clear();
//	std::vector<int> finalSegSeedIndexes;
//
//	
//	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
//	{
//		segmentsSize.insert({ bigSegIt->first, static_cast<double>(bigSegIt->second.piece_vertices_index_.size()) });
//		segmentsAvgNormal.insert(
//			{ bigSegIt->first, fragment.localAvgNormal(bigSegIt->second.piece_vertices_index_).normalized() }
//		);
//		segmentsAvgCurvedness.insert(
//			{ bigSegIt->first,fragment.localAvgCurvedness(bigSegIt->second.piece_vertices_index_) }
//		);
//
//		Eigen::Vector3d sums = Eigen::Vector3d::Zero();
//
//		for (int iVert : bigSegIt->second.piece_vertices_index_)
//		{
//			for (int j = 0; j < 3; j++)
//			{
//				sums[j] += fragment.m_Vertices.coeff(iVert, j);
//			}
//
//		}
//
//		segmentsCenterMass.insert({ bigSegIt->first, (sums / bigSegIt->second.piece_vertices_index_.size()).normalized()});
//	}
//
//	int iIntactSeg = std::min_element(
//						segmentsAvgCurvedness.begin(), segmentsAvgCurvedness.end(),
//						[](const auto& l, const auto& r) {return l.second < r.second; })->first;
//	int iBiggestSeg = std::max_element(
//						segmentsSize.begin(), segmentsSize.end(),
//						[](const auto& l, const auto& r) {return l.second < r.second; })->first;
//
//	if (iIntactSeg != iBiggestSeg)
//	{
//		std::cout << "Note: the intact surface is not the biggest: check the result of the segmentation" << std::endl;
//	}
//	
//	//finalSegSeedIndexes.push_back(iIntactSeg);
//
//	/*
//		Find the seed of the opposite segment
//	*/
//
//	double minSim = 2;
//	int iOppositeSeg = -1;
//	for (auto bigSegIt = segmentsAvgNormal.begin(); bigSegIt != segmentsAvgNormal.end(); bigSegIt++)
//	{
//		double sim = bigSegIt->second.dot(segmentsAvgNormal.at(iIntactSeg));
//
//		if (sim < minSim)
//		{
//			iOppositeSeg = bigSegIt->first;
//			minSim = sim;
//		}
//	}
//
//	finalSegSeedIndexes.push_back(iOppositeSeg);
//
//
//	/*
//		Find the side walls and opposite segment seeds
//	*/
//	std::vector<int> sidewallsSegIndexes;
//	std::vector<double> wallsSegSimToIntact;
//	double eplisionOrthErr = 0.1; //0.1; // it depends on the size of the segments
//	
//	for (auto bigSegIt = segmentsAvgNormal.begin(); bigSegIt != segmentsAvgNormal.end(); bigSegIt++)
//	{
//
//		if (bigSegIt->first == iIntactSeg || bigSegIt->first == iOppositeSeg)
//		{
//			continue;
//		}
//
//		double simToIntact = abs(bigSegIt->second.dot(segmentsAvgNormal.at(iIntactSeg)));
//
//		if (simToIntact > eplisionOrthErr)
//		{
//			continue;
//		}
//
//		sidewallsSegIndexes.push_back(bigSegIt->first);
//
//
//
//	}
//
//	
//
//	std::vector< Eigen::RowVector3d> dirs = {
//		{0,0,1},
//		{0,0,-1},
//		{0,1,0},
//		{0,-1,0},
//		{1,0,0},
//		{-1,0,0}
//	};
//
//	std::map<int, std::pair<int,double>> dir2walls;
//
//	for (int iWallSeg: sidewallsSegIndexes)
//	{
//		
//		auto wallOrientation = segmentsCenterMass.at(iWallSeg) - segmentsCenterMass.at(iIntactSeg);
//		int iFitDir = -1;
//		double maxSim = -2;
//		double sim = -3;
//
//		for (int iDir = 0; iDir < dirs.size();iDir++)
//		{
//			sim = dirs[iDir].dot(wallOrientation);
//			if (sim > maxSim)
//			{
//				maxSim = sim;
//				iFitDir = iDir;
//			}
//
//		}
//		
//		if (dir2walls.count(iFitDir)==0)
//		{
//			dir2walls[iFitDir] = { iWallSeg,maxSim };
//		}
//		else
//		{
//			if (dir2walls.at(iFitDir).second < maxSim)
//			{
//
//				dir2walls[iFitDir].first =  iWallSeg;
//				dir2walls[iFitDir].second =  maxSim;
//			}
//		}
//		
//	}
//
//	for (auto &itSeg : dir2walls)
//	{
//		finalSegSeedIndexes.push_back(itSeg.second.first);
//	}
//
//	Eigen::MatrixXd finalSegSeedColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
//	colorIt = meshColors.begin();
//
//	for (auto it = finalSegSeedIndexes.begin(); it != finalSegSeedIndexes.end(); it++)
//	{
//		
//		int iSeg =*it;
//			std::vector<int> regionVerticesIndx = bigSegments.at(iSeg).piece_vertices_index_;
//			for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//			{
//				finalSegSeedColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//			}
//		
//		++colorIt;
//	}
//
//	
//
//	eplisionOrthErr = 0.1;
//	std::map<int, std::pair<int,double>> seg2Seeds;
//
//
//	while (bigSegments.size() > 6) 
//	{
//
//
//		for (int iSegFinal : finalSegSeedIndexes)
//		{
//			auto& iVertsAtBoundary = oRegionOutsideBoundaryVerticesList[iSegFinal];
//
//			for (int iVertBoundary : iVertsAtBoundary)
//			{
//				int iNeighboorSeg = vertIndex2SegIndex[iVertBoundary];
//
//				// Don't mix between seeds
//				auto& itNeigh = std::find(finalSegSeedIndexes.begin(), finalSegSeedIndexes.end(), iNeighboorSeg);
//
//				if (itNeigh!=finalSegSeedIndexes.end())
//				{
//					continue;
//				}
//
//				if (iNeighboorSeg != iSegFinal && iNeighboorSeg != iIntactSeg)
//				{
//					double currSimilariry = segmentsAvgNormal[iNeighboorSeg].dot(segmentsAvgNormal[iSegFinal]);
//
//					if (1 - currSimilariry < eplisionOrthErr)
//					{
//						if (seg2Seeds.count(iNeighboorSeg) == 0)
//						{
//							seg2Seeds[iNeighboorSeg] = { iSegFinal,currSimilariry };
//						}
//						else
//						{
//							if (currSimilariry > seg2Seeds[iNeighboorSeg].second)
//							{
//								seg2Seeds[iNeighboorSeg].first = iSegFinal;
//								seg2Seeds[iNeighboorSeg].second = currSimilariry;
//							}
//						}
//					}
//
//				}
//			}
//
//		}
//
//		for (auto& itSeg = seg2Seeds.begin(); itSeg != seg2Seeds.end(); itSeg++)
//		{
//			int iSegDst = itSeg->second.first;
//			int iSegSrc = itSeg->first;
//			for (int iVert : bigSegments.at(iSegSrc).piece_vertices_index_)
//			{
//				vertIndex2SegIndex[iVert] = iSegDst;
//			}
//
//
//			bigSegments.at(iSegDst).piece_vertices_index_.insert(
//				bigSegments.at(iSegDst).piece_vertices_index_.end(),
//				bigSegments.at(iSegSrc).piece_vertices_index_.begin(),
//				bigSegments.at(iSegSrc).piece_vertices_index_.end()
//			);
//
//			oRegionOutsideBoundaryVerticesList[iSegDst].insert(
//				oRegionOutsideBoundaryVerticesList[iSegDst].end(),
//				oRegionOutsideBoundaryVerticesList[iSegSrc].begin(),
//				oRegionOutsideBoundaryVerticesList[iSegSrc].end()
//			);
//
//
//			segmentsSize[iSegDst] = static_cast<double>(bigSegments[iSegDst].piece_vertices_index_.size());
//			segmentsAvgNormal[iSegDst] = fragment.localAvgNormal(bigSegments[iSegDst].piece_vertices_index_).normalized();
//			bigSegments.erase(iSegSrc);
//			segmentsAvgNormal.erase(iSegSrc);
//			segmentsSize.erase(iSegSrc);
//		}
//
//		if (seg2Seeds.size()==0)
//		{
//			eplisionOrthErr += 0.02;
//
//		}
//		seg2Seeds.clear();
//	}
//
//	Eigen::MatrixXd finalSegSeedColorsTmp = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
//	colorIt = meshColors.begin();
//
//	finalSegSeedIndexes.push_back(iIntactSeg);
//
//	for (auto iSeg: finalSegSeedIndexes)
//	{
//
//
//		std::vector<int> regionVerticesIndx = bigSegments.at(iSeg).piece_vertices_index_;
//		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
//		{
//			finalSegSeedColorsTmp.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
//		}
//
//		++colorIt;
//	}
//	
//
//
//	visualizer.m_Viewer.callback_key_down =
//		[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) ->bool
//	{
//
//		switch (key) {
//		case '1': // for debug
//			visualizer.m_Viewer.data().set_colors(fragment.m_Colors);//meshColors[0]
//			std::cout << "Pressed 1" << std::endl;
//			break;
//		case '2': // for debug
//			// before merging
//			visualizer.m_Viewer.data().set_colors(rawSegmentColors);
//			std::cout << "Pressed 2, present the small and big regions before merge" << std::endl;
//			//visualizer.writeOFF("beforeMerge.off", fragment.m_Vertices, fragment.m_Faces, segmentColors.block(0,0, segmentColors.rows(),3));
//			break;
//		case '3': // for debug
//			// Only big segments
//			visualizer.m_Viewer.data().set_colors(bigSegmentColors);
//			std::cout << "Pressed 3, present only the big segments" << std::endl;
//			break;
//		case '4': // for debug
//			visualizer.m_Viewer.data().set_colors(segmentColorsAfterMerge);
//			std::cout << "Pressed 4, present the mesh after merging" << std::endl;
//			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
//			break;
//
//
//		case '5': // for debug
//			visualizer.m_Viewer.data().set_colors(segmentBoundaryAfterMerge);
//			std::cout << "Pressed 5, present the mesh boundary after merging" << std::endl;
//			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
//			break;
//
//		case '6': // for debug
//			visualizer.m_Viewer.data().set_colors(finalSegSeedColors);
//			std::cout << "Pressed 6, show the seeds for the final segmentse" << std::endl;
//			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
//			break;
//		
//		case '7': // for debug
//			visualizer.m_Viewer.data().set_colors(finalSegSeedColorsTmp);
//			std::cout << "Pressed 7,finalSegSeedColorsTmp" << std::endl;
//			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
//			break;
//		}
//		
//		return false;
//	};
//	visualizer.launch();
//
//	//std::vector<Segment> segments;
//	//std::vector<Segment> segments; // make it call by reference
//	//for (int iRegion = 0; iRegion < oRegionsList.size(); iRegion++)
//	//{
//	//	Segment segment(oRegionsList[iRegion]);
//	//	segment.loadBasicData(fragment.m_VerticesAdjacentFacesList,
//	//		fragment.m_Faces, fragment.m_Vertices,
//	//		fragment.m_Faces2TextureCoordinates, fragment.m_Faces2Normals);
//	//	segments.push_back(segment);
//	//}
//}