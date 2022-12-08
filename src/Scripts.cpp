#include "Scripts.h"
#include "Visualization.h"


void segment_intact_surface(std::vector<std::string> all_args)
{
	std::string fragmentPath = all_args[0];
	std::string outFileName = all_args[1];


	ObjFragment fragment = ObjFragment(fragmentPath);
	std::cout << "Loading data" << std::endl;
	fragment.load();

	std::cout << "Start to segment" << std::endl;
	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;

	bool isSegmented = false;

	double fracture = 0.65;
	std::vector<Segment> segments;
	int nTrials = 1;

	while (!isSegmented)
	{
		double simThresh = fragment.getSimilarThreshByPos(fracture);
		std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
		fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);
		fragment.filterSmallRegions(segments, oRegionsList);
		std::cout << "In trial: " << nTrials << " Found " << segments.size() << " Segments" << std::endl;

		if (segments.size() == 0)
		{
			fracture = fracture + 0.05;
		}
		else {
			if (segments.size() == 1)
			{
				fracture = fracture - 0.12;
			}
			else
			{
				isSegmented = true;
			}

		}

		nTrials += 1;

		if (fracture > 1)
		{
			std::cout << "fracrure value is " << fracture << " and it should be between 0 to 1" << std::endl;
			std::cout << "Please rerun mannually..exiting" << std::endl;
			std::exit(1);
		}

	}

	int intactIndex = fragment.findIntactSegmentIndex(segments);

	Segment intactSurface(segments[intactIndex]);
	intactSurface.loadBasicData(fragment.m_VerticesAdjacentFacesList,
								fragment.m_Faces, fragment.m_Vertices,
								fragment.m_Faces2TextureCoordinates, fragment.m_Faces2Normals);

	std::string fragFolderPath = fragmentPath.substr(0, fragmentPath.find_last_of("\\/"));
	fragment.saveAsObj(fragFolderPath + "\\" + outFileName, intactSurface);
	std::cout << "Write successfully the output to path " << fragFolderPath << std::endl;


	// For debug:
	
}

void segment(std::vector<std::string> all_args)
{
	std::string fragmentPath = all_args[0];
	std::string outFileName = all_args[1];
	ObjFragment fragment = ObjFragment(fragmentPath);
	std::cout << "Loading data" << std::endl;
	fragment.load();

	std::cout << "Start to segment" << std::endl;
	double fracture = 0.5; //0.45; //0.5
	double simThresh = fragment.getSimilarThreshByPos(fracture);
	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);


	/*
	
		Start mergnig algo - you can put it in other function when refactor
	*/

	/*
	* init data structures
	*/
	double minSegPercSize = 0.000125;  //0.000125;//0.005 // This need to get as a parameter to the script
	std::map<int, Segment> smallSegments;
	std::map<int, Segment> bigSegments;
	std::map<int, double> segmentsAvgCurvedness; // We avg? maybe std is a parameter we should consider
	
	std::map<int, int> vertIndex2SegIndex;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	std::map<int, int> freeVertIndexes;

	for (int iSeg = 0; iSeg < oRegionsList.size(); iSeg++)
	{
		auto segVertsIndexes = oRegionsList[iSeg];
		Segment currSeg(segVertsIndexes);
		segmentsAvgCurvedness.insert({ iSeg,fragment.localAvgCurvedness(segVertsIndexes) });
		

		double segmentSize = static_cast<double>(currSeg.piece_vertices_index_.size());
		if (segmentSize/fragmentSize < minSegPercSize){
			smallSegments.insert({ iSeg, currSeg });
			

			for (int iVert : segVertsIndexes)
			{
				freeVertIndexes.insert({ iVert,-1 });
			}
		}
		else
		{
			bigSegments.insert({ iSeg, currSeg });

		}

		for (int i = 0; i < segVertsIndexes.size(); i++)
		{
			vertIndex2SegIndex.insert({ segVertsIndexes[i],iSeg });
		}

	}

	/*
		For debug
	*/
	int ixMaxSize_debug = -1;
	double maxSizeSegment_debug = 0;
	int k = 0;
	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	{
		double segmentSize = static_cast<double>(bigSegIt->second.piece_vertices_index_.size());
		if (segmentSize > maxSizeSegment_debug)
		{
			maxSizeSegment_debug = segmentSize;
			ixMaxSize_debug = k;
		}
		++k;
	}

	/*
		Coloring
	*/
	Visualizer visualizer;
	std::map<int, Eigen::RowVector3d> meshColors;
	visualizer.generateRandomColors(meshColors, oRegionsList.size());
	visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, meshColors[0]);

	Eigen::MatrixXd bigSegmentColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);

	auto colorIt = meshColors.begin();
	//for (auto bigSegIt = std::next(bigSegments.begin(), ixMaxSize_debug); bigSegIt != std::next(bigSegments.begin(), ixMaxSize_debug + 1); bigSegIt++)
	for (auto bigSegIt = bigSegments.begin(); bigSegIt !=bigSegments.end(); bigSegIt++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[bigSegIt->first];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{

			bigSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		++colorIt;
	}

	Eigen::MatrixXd rawSegmentColors;
	rawSegmentColors.resize(fragment.m_Vertices.rows(), 4);
	//segmentBigSmallColors.resize(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();
	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[bigSegIt->first];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			//bigSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
			rawSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		
			/*visualizer.m_Viewer.data().add_label(fragment.m_Vertices.row(regionVerticesIndx[iVert]),
													std::to_string(bigSegIt->first));*/
		}

		colorIt++;
	}

	for (auto smallSegIt = smallSegments.begin(); smallSegIt != smallSegments.end(); smallSegIt++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[smallSegIt->first];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			rawSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		colorIt++;
	}


	std::set<int> ixCurrSegBigNeigh;
	std::set<int> ixCurrSegSmallNeigh;
	std::map<int, int> segSrc2segDst; // The src seg that will be merge to dst seg
	std::map<int, std::vector<int>> segSrc2MultiDst;

	std::cout << "Start the merge process" << std::endl;
	int nLastFreeDebug = -1;
	int debugIter = 0;
	std::cout << "WARNING: you run only one iteration" << std::endl;
	while (!smallSegments.empty()) //debugIter++ < 1()
	{
		debugIter++;
		if (nLastFreeDebug == smallSegments.size())
		{
			std::cout << "WARNING: merge did not converged after " << debugIter << " iterations" << std::endl;
			break;
		}
		nLastFreeDebug = smallSegments.size();

		//std::cout << "WARNING: merge only for a single big segment" << std::endl;
		//for (auto bigSegIt = std::next(bigSegments.begin(), ixMaxSize_debug); bigSegIt != std::next(bigSegments.begin(), ixMaxSize_debug+1); bigSegIt++)//bigSegments.end()
		for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)//bigSegments.end()
		{
			auto& iVertsAtBoundary = oRegionOutsideBoundaryVerticesList[bigSegIt->first];
			
			/*
				Compute the neighboors segments
			*/
			for (int iVertBoundary : iVertsAtBoundary)
			{
				int iNeighboorSeg = vertIndex2SegIndex[iVertBoundary];

				// Identify wheter its a big or small segment
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
				//segSrc2segDst.insert({ *itSmallNeigh,bigSegIt->first});
				//segSrc2MultiDst[*itSmallNeigh].push_back
			}
			

			ixCurrSegBigNeigh.clear();
			ixCurrSegSmallNeigh.clear();
		}

		//std::cout << "before applying the merging" << std::endl;
		/*
			Merging
		*/
		for (std::map<int,int>::iterator itSegSrc2Dst = segSrc2segDst.begin(); itSegSrc2Dst !=segSrc2segDst.end();++itSegSrc2Dst)
		{
			int iSrcSeg = itSegSrc2Dst->first;
			int iDstSeg = itSegSrc2Dst->second;

			
			if (smallSegments.count(iSrcSeg) > 0)
			{


				for (int i = 0; i < smallSegments.at(iSrcSeg).piece_vertices_index_.size(); i++)
				{
					vertIndex2SegIndex[smallSegments.at(iSrcSeg).piece_vertices_index_[i]] = iDstSeg;
				}

				if (bigSegments.count(iDstSeg) > 0)
				{
					bigSegments.at(iDstSeg).piece_vertices_index_.insert(
						bigSegments.at(iDstSeg).piece_vertices_index_.end(),
						smallSegments.at(iSrcSeg).piece_vertices_index_.begin(),
						smallSegments.at(iSrcSeg).piece_vertices_index_.end()
					);

					oRegionOutsideBoundaryVerticesList[iDstSeg].insert(
						oRegionOutsideBoundaryVerticesList[iDstSeg].end(),
						oRegionOutsideBoundaryVerticesList[iSrcSeg].begin(),
						oRegionOutsideBoundaryVerticesList[iSrcSeg].end()
					);
				}
				else {
					std::cout << "You should not be here";
				}

			}
			
			
		}

		for (std::map<int, int>::iterator itSrc2Dst = segSrc2segDst.begin(); itSrc2Dst != segSrc2segDst.end(); ++itSrc2Dst)
		{
			/*auto currSeg = &smallSegments.at(itSrc2Dst->first);
			double segmentSize = static_cast<double>(currSeg->piece_vertices_index_.size());
			if (segmentSize / fragmentSize > minSegPercSize) {
				bigSegments.insert({ itSrc2Dst->first ,*currSeg });
			}*/
			smallSegments.erase(itSrc2Dst->first);
		}
		segSrc2segDst.clear();
	}
	
	/*
		Coloring
	*/

	Eigen::MatrixXd segmentColorsAfterMerge = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();

	//std::cout << "WARNING: You present the merge only for a single big segment" << std::endl;
	//for (auto bigSegIt = std::next(bigSegments.begin(), ixMaxSize_debug); bigSegIt != std::next(bigSegments.begin(), ixMaxSize_debug + 1); bigSegIt++)
	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	{
		std::vector<int> regionVerticesIndx = bigSegIt->second.piece_vertices_index_;
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			segmentColorsAfterMerge.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		++colorIt;
	}

	/*
		Coloring
	*/
	Eigen::MatrixXd segmentBoundaryAfterMerge = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();

	for (int iRegion = 0; iRegion < oRegionOutsideBoundaryVerticesList.size(); iRegion++)
	{
		if (bigSegments.count(iRegion) > 0)
		{

			std::vector<int> regionVerticesIndx = oRegionOutsideBoundaryVerticesList[iRegion];
			for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
			{
				segmentBoundaryAfterMerge.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
			}
			++colorIt;
		}
	}

	/*
		Merge the big segments by normals 
	*/

	if (bigSegments.size() < 6)
	{
		std::cout << "ERROR: Expect that the final number of big segment to be greater than 6 (Prior Knoledge)" << std::endl;
		exit(1);
	}

	std::map<int, Eigen::Vector3d> segmentsAvgNormal;
	std::map<int, double> segmentsSize;
	segmentsAvgCurvedness.clear();
	std::vector<int> finalSegSeedIndexes;

	
	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	{
		segmentsSize.insert({ bigSegIt->first, static_cast<double>(bigSegIt->second.piece_vertices_index_.size()) });
		segmentsAvgNormal.insert(
			{ bigSegIt->first, fragment.localAvgNormal(bigSegIt->second.piece_vertices_index_).normalized() }
		);
		segmentsAvgCurvedness.insert(
			{ bigSegIt->first,fragment.localAvgCurvedness(bigSegIt->second.piece_vertices_index_) }
		);
	}

	int iIntactSeg = std::min_element(
						segmentsAvgCurvedness.begin(), segmentsAvgCurvedness.end(),
						[](const auto& l, const auto& r) {return l.second < r.second; })->first;
	int iBiggestSeg = std::max_element(
						segmentsSize.begin(), segmentsSize.end(),
						[](const auto& l, const auto& r) {return l.second < r.second; })->first;

	if (iIntactSeg != iBiggestSeg)
	{
		std::cout << "Note: the intact surface is not the biggest: check the result of the segmentation" << std::endl;
	}
	
	finalSegSeedIndexes.push_back(iIntactSeg);

	/*
		Find the seed of the opposite segment
	*/

	double minSim = 2;
	int iOppositeSeg = -1;
	for (auto bigSegIt = segmentsAvgNormal.begin(); bigSegIt != segmentsAvgNormal.end(); bigSegIt++)
	{
		double sim = bigSegIt->second.dot(segmentsAvgNormal.at(iIntactSeg));

		if (sim < minSim)
		{
			iOppositeSeg = bigSegIt->first;
			minSim = sim;
		}
	}

	finalSegSeedIndexes.push_back(iOppositeSeg);


	/*
		Find the side walls and opposite segment seeds
	*/
	std::vector<int> sidewallsSegIndexes;
	sidewallsSegIndexes.push_back(iIntactSeg);
	double eplisionOrthErr = 0.00001;
	for (auto bigSegIt = segmentsAvgNormal.begin(); bigSegIt != segmentsAvgNormal.end(); bigSegIt++)
	{
		//int nOrthFracSeg = 0;
		//int iSegAlign = -1;
		//double sim = -2;

		//for (int i=0; i<sidewallsSegIndexes.size();++i) 
		//{
		//	int iFracSeg = sidewallsSegIndexes[i];
		//	sim = bigSegIt->second.dot(segmentsAvgNormal.at(iFracSeg));

		//	// If it is orthogonal
		//	if (sim < eplisionOrthErr)
		//	{
		//		nOrthFracSeg++;
		//	}

		//	// If it is the same direction to one of the vector orthogonal to the intact
		//	if (1-sim < eplisionOrthErr)
		//	{
		//		iSegAlign = iFracSeg;
		//	}
		//}

		//if (nOrthFracSeg < sidewallsSegIndexes.size() - 1 )
		//{
		//	if (iSegAlign == -1)
		//	{
		//		sidewallsSegIndexes.push_back(bigSegIt->first);
		//	}
		//	else
		//	{

		//	}
		//}

	}

	Eigen::MatrixXd finalSegSeedColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();

	for (int iSeg: finalSegSeedIndexes)
	{
		std::vector<int> regionVerticesIndx = bigSegments.at(iSeg).piece_vertices_index_;
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			finalSegSeedColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		++colorIt;
	}



	//int iMostSimSegNeigh=-1;
	//double maxSim = -2;
	//int nDebug = 0;

	//while (bigSegments.size() > 6) // nDebug++ < 1
	//{
	//	int iMinSegment = std::min_element(segmentsSize.begin(), segmentsSize.end(), [](const auto& l, const auto& r) {return l.second < r.second; })->first;
	//	auto& iVertsAtBoundary = oRegionOutsideBoundaryVerticesList[iMinSegment];

	//	iMostSimSegNeigh = -1;
	//	maxSim = -2;

	//	/*
	//		Compute the neighboors segments
	//	*/
	//	for (int iVertBoundary : iVertsAtBoundary)
	//	{
	//		int iNeighboorSeg = vertIndex2SegIndex[iVertBoundary];

	//		if (iNeighboorSeg != iMinSegment)
	//		{
	//			double currSimilariry = segmentsAvgNormal[iNeighboorSeg].dot(segmentsAvgNormal[iMinSegment]);

	//			if (currSimilariry > maxSim)
	//			{
	//				iMostSimSegNeigh = iNeighboorSeg;
	//				maxSim = currSimilariry;
	//			}
	//		}
	//	}

	//	if (iMostSimSegNeigh==-1)
	//	{
	//		std::cout << "In the second merge: we must have fit neighboors" << std::endl;
	//		std::exit(1);
	//	}

	//	/*for (int i = 0; i < bigSegments.at(iMinSegment).piece_vertices_index_.size(); i++)
	//	{
	//		vertIndex2SegIndex[bigSegments.at(iMinSegment).piece_vertices_index_[i]] = iMostSimSegNeigh;
	//	}*/
	//	for (int iVert : bigSegments.at(iMinSegment).piece_vertices_index_)
	//	{
	//		vertIndex2SegIndex[iVert] = iMostSimSegNeigh;
	//	}

	//		
	//	bigSegments.at(iMostSimSegNeigh).piece_vertices_index_.insert(
	//		bigSegments.at(iMostSimSegNeigh).piece_vertices_index_.end(),
	//		bigSegments.at(iMinSegment).piece_vertices_index_.begin(),
	//		bigSegments.at(iMinSegment).piece_vertices_index_.end()
	//	);

	//	oRegionOutsideBoundaryVerticesList[iMostSimSegNeigh].insert(
	//		oRegionOutsideBoundaryVerticesList[iMostSimSegNeigh].end(),
	//		oRegionOutsideBoundaryVerticesList[iMinSegment].begin(),
	//		oRegionOutsideBoundaryVerticesList[iMinSegment].end()
	//	);


	//	segmentsSize[iMostSimSegNeigh]= static_cast<double>(bigSegments[iMostSimSegNeigh].piece_vertices_index_.size());
	//	segmentsAvgNormal[iMostSimSegNeigh] = fragment.localAvgNormal(bigSegments[iMostSimSegNeigh].piece_vertices_index_).normalized();
	//	bigSegments.erase(iMinSegment);
	//	segmentsAvgNormal.erase(iMinSegment);
	//	segmentsSize.erase(iMinSegment);
	//}

	///*
	//	Coloring

	//*/

	//Eigen::MatrixXd segmentAfterSecondMerge = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	//colorIt = meshColors.begin();

	//for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	//{
	//	//auto bigSegIt = bigSegments.find(iMostSimSegNeigh);
	//	std::vector<int> regionVerticesIndx = bigSegIt->second.piece_vertices_index_;
	//	for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
	//	{
	//		segmentAfterSecondMerge.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
	//	}
	//	++colorIt;
	//}

	// meshColors[0] //segmentColors
	visualizer.m_Viewer.callback_key_down =
		[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) ->bool
	{

		switch (key) {
		case '1': // for debug
			visualizer.m_Viewer.data().set_colors(meshColors[0]);
			std::cout << "Pressed 1" << std::endl;
			break;
		case '2': // for debug
			// before merging
			visualizer.m_Viewer.data().set_colors(rawSegmentColors);
			std::cout << "Pressed 2, present the small and big regions before merge" << std::endl;
			//visualizer.writeOFF("beforeMerge.off", fragment.m_Vertices, fragment.m_Faces, segmentColors.block(0,0, segmentColors.rows(),3));
			break;
		case '3': // for debug
			// Only big segments
			visualizer.m_Viewer.data().set_colors(bigSegmentColors);
			std::cout << "Pressed 3, present only the big segments" << std::endl;
			break;
		case '4': // for debug
			visualizer.m_Viewer.data().set_colors(segmentColorsAfterMerge);
			std::cout << "Pressed 4, present the mesh after merging" << std::endl;
			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
			break;


		case '5': // for debug
			visualizer.m_Viewer.data().set_colors(segmentBoundaryAfterMerge);
			std::cout << "Pressed 5, present the mesh boundary after merging" << std::endl;
			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
			break;

		case '6': // for debug
			visualizer.m_Viewer.data().set_colors(finalSegSeedColors);
			std::cout << "Pressed 6, show the seeds for the final segmentse" << std::endl;
			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
			break;
		
		//case '7': // for debug
		//	visualizer.m_Viewer.data().set_colors(segmentAfterSecondMerge);
		//	std::cout << "Pressed 6, present the mesh after second merging" << std::endl;
		//	//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
		//	break;
		}
		
		return false;
	};
	visualizer.launch();

	//std::vector<Segment> segments;
	//std::vector<Segment> segments; // make it call by reference
	//for (int iRegion = 0; iRegion < oRegionsList.size(); iRegion++)
	//{
	//	Segment segment(oRegionsList[iRegion]);
	//	segment.loadBasicData(fragment.m_VerticesAdjacentFacesList,
	//		fragment.m_Faces, fragment.m_Vertices,
	//		fragment.m_Faces2TextureCoordinates, fragment.m_Faces2Normals);
	//	segments.push_back(segment);
	//}
}