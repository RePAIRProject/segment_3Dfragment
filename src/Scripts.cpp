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
	double fracture = 0.5;
	double simThresh = fragment.getSimilarThreshByPos(fracture);
	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);

	Visualizer visualizer;
	std::map<int, Eigen::RowVector3d> meshColors;
	visualizer.generateRandomColors(meshColors, oRegionsList.size());

	Eigen::MatrixXd segmentColors; 
	segmentColors.resize(fragment.m_Vertices.rows(),4);
	auto colorIt = meshColors.begin();

	for (int iRegion = 0; iRegion < oRegionsList.size(); iRegion++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[iRegion];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			segmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		++colorIt;
	}



	/*
	
		Start mergnig algo - you can put it in other function when refactor
	*/

	/*
	* init data structures
	*/
	double minSegPercSize = 0.005; // This need to get as a parameter to the script
	std::map<int, Segment> smallSegments;
	std::map<int, Segment> bigSegments;
	std::map<int, double> segmentsAvgCurvedness; // We avg? maybe std is a parameter we should consider
	std::map<int, Eigen::Vector3d> segmentsAvgNormal;
	std::map<int, int> vertIndex2SegIndex;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	for (int iSeg = 0; iSeg < oRegionsList.size(); iSeg++)
	{
		auto segVertsIndexes = oRegionsList[iSeg];
		Segment currSeg(segVertsIndexes);
		segmentsAvgCurvedness.insert({ iSeg,fragment.localCurvedness(segVertsIndexes) });
		segmentsAvgNormal.insert({ iSeg, fragment.localAvgNormal(segVertsIndexes).normalized()});

		double segmentSize = static_cast<double>(currSeg.piece_vertices_index_.size());
		if (segmentSize/fragmentSize < minSegPercSize){
			smallSegments.insert({ iSeg, currSeg });
		}
		else
		{
			bigSegments.insert({ iSeg, currSeg });
		}

		for (int iVert = 0; iVert < segVertsIndexes.size(); iVert++)
		{
			vertIndex2SegIndex.insert({ segVertsIndexes[iVert],iSeg });
		}

	}


	Eigen::MatrixXd segmentBigSmallColors;
	segmentBigSmallColors.resize(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();

	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[bigSegIt->first];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			segmentBigSmallColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}

		if (colorIt == std::next(meshColors.begin(), 1))
		{
			colorIt = meshColors.begin();
		}
		else 
		{
			++colorIt;
		}
	}

	colorIt = meshColors.begin();

	for (auto smallSegIt = smallSegments.begin(); smallSegIt != smallSegments.end(); smallSegIt++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[smallSegIt->first];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			segmentBigSmallColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		++colorIt;
	}

	/*
		Find the neighboorhod segments of the current segment
	*/
	//std::map<int, std::set<int>> seg2BigNeigh;
	//std::map<int, std::set<int>> seg2SmallNeigh;
	//std::map<int, std::set<int>>* smallOrBigMap;
	/*std::set<int> ixCurrSegBigNeigh;
	std::set<int> ixCurrSegSmallNeigh;


	for (auto smallSegIt = smallSegments.begin(); smallSegIt != smallSegments.end(); smallSegIt++)
	{
		std::set<int> emptySet;
		seg2BigNeigh.insert({ smallSegIt->first,emptySet });
		seg2SmallNeigh.insert({ smallSegIt->first,emptySet });
	}*/

	std::set<int> ixCurrSegBigNeigh;
	std::set<int> ixCurrSegSmallNeigh;
	std::map<int, int> segSrc2segDst; // The src seg that will be merge to dst seg

	while (!smallSegments.empty())
	{
		for (auto smallSegIt = smallSegments.begin(); smallSegIt != smallSegments.end(); smallSegIt++)
		{
			auto& iVertsAtBoundary = oRegionOutsideBoundaryVerticesList[smallSegIt->first];
			
			/*
				Compute the neighboors segments
			*/
			for (int iVertBoundary = 0; iVertBoundary < iVertsAtBoundary.size(); iVertBoundary++)
			{
				int iNeighboorSeg = vertIndex2SegIndex[iVertBoundary];

				// Identify wheter its a big or small segment
				if (smallSegments.count(iNeighboorSeg) > 0)
				{
					ixCurrSegSmallNeigh.insert(iNeighboorSeg);
				}
				else
				{
					ixCurrSegBigNeigh.insert(iNeighboorSeg);
				}
			}

			/*
				Find the most similar in the big segments nieghbors
			*/
			bool isMerged = false;

			if (!ixCurrSegBigNeigh.empty())
			{

				/*int iBigSegMostSim = -1;
				double minCurvDiff = 99999;

				for (std::set<int>::iterator itBigNeigh = ixCurrSegBigNeigh.begin(); itBigNeigh != ixCurrSegBigNeigh.end(); itBigNeigh++)
				{

					double currDiff = abs(segmentsAvgCurvedness.at(*itBigNeigh) - segmentsAvgCurvedness.at(smallSegIt->first));
					if (currDiff < minCurvDiff)
					{
						iBigSegMostSim = *itBigNeigh;
						minCurvDiff = currDiff;
					}
				}*/

				int iBigSegMostSim = -1;
				double maxInnerProd = -2; // The vectors should be normalized

				for (std::set<int>::iterator itBigNeigh = ixCurrSegBigNeigh.begin(); itBigNeigh != ixCurrSegBigNeigh.end(); itBigNeigh++)
				{

					double currInnerProd = segmentsAvgNormal.at(*itBigNeigh).dot(segmentsAvgNormal.at(smallSegIt->first));
					if (currInnerProd > maxInnerProd)
					{
						iBigSegMostSim = *itBigNeigh;
						maxInnerProd = currInnerProd;
					}
				}

				if (iBigSegMostSim != -1)
				{
					segSrc2segDst.insert({ smallSegIt->first,iBigSegMostSim });
				}
			}
			
			ixCurrSegBigNeigh.clear();
			ixCurrSegSmallNeigh.clear();
		}

		/*
			Merging
		*/
		for (std::map<int,int>::iterator itsmall2Big = segSrc2segDst.begin(); itsmall2Big !=segSrc2segDst.end();++itsmall2Big)
		{
			int iSmallSeg = itsmall2Big->first;
			int iBigSeg = itsmall2Big->second;

			for (int iVert = 0; iVert < smallSegments.at(iSmallSeg).piece_vertices_index_.size(); iVert++)
			{
				vertIndex2SegIndex[iVert] = iBigSeg;
			}

			bigSegments.at(iBigSeg).piece_vertices_index_.insert(
						bigSegments.at(iBigSeg).piece_vertices_index_.end(), 
						smallSegments.at(iSmallSeg).piece_vertices_index_.begin(),
						smallSegments.at(iSmallSeg).piece_vertices_index_.end()
						);

			oRegionOutsideBoundaryVerticesList[iBigSeg].insert(
				oRegionOutsideBoundaryVerticesList[iBigSeg].end(),
				oRegionOutsideBoundaryVerticesList[iSmallSeg].begin(),
				oRegionOutsideBoundaryVerticesList[iSmallSeg].end()
			);

		}

		for (std::map<int, int>::iterator itsmall2Big = segSrc2segDst.begin(); itsmall2Big != segSrc2segDst.end(); ++itsmall2Big)
		{
			smallSegments.erase(itsmall2Big->first);
		}
		segSrc2segDst.clear();
	}
	
	Eigen::MatrixXd segmentColorsAfterMerge;
	segmentColorsAfterMerge.resize(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();

	for (auto it = bigSegments.begin(); it != bigSegments.end(); ++it)
	{
		std::vector<int> regionVerticesIndx = it->second.piece_vertices_index_;
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			segmentColorsAfterMerge.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		++colorIt;
	}


	visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, segmentColors);// meshColors[0] //segmentColors
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
			visualizer.m_Viewer.data().set_colors(segmentColors);
			std::cout << "Pressed 2" << std::endl;
			//visualizer.writeOFF("beforeMerge.off", fragment.m_Vertices, fragment.m_Faces, segmentColors.block(0,0, segmentColors.rows(),3));
			break;
		case '3': // for debug
			// before merging
			visualizer.m_Viewer.data().set_colors(segmentBigSmallColors);
			std::cout << "Pressed 3" << std::endl;
			break;
		case '4': // for debug
			// before merging
			visualizer.m_Viewer.data().set_colors(segmentColorsAfterMerge);
			std::cout << "Pressed 4" << std::endl;
			//visualizer.writeOFF("byNormals.off", fragment.m_Vertices, fragment.m_Faces, segmentColorsAfterMerge.block(0, 0, segmentColors.rows(), 3));
			break;
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