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
	double fracture = 0.45; //0.5
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
	double minSegPercSize = 0.000125;//0.005 // This need to get as a parameter to the script
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

	/*
		Coloring
	*/

	Eigen::MatrixXd bigSegmentColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	//segmentBigSmallColors.resize(fragment.m_Vertices.rows(), 4);
	colorIt = meshColors.begin();

	for (auto bigSegIt = bigSegments.begin(); bigSegIt != bigSegments.end(); bigSegIt++)
	{
		std::vector<int> regionVerticesIndx = oRegionsList[bigSegIt->first];
		for (int iVert = 0; iVert < regionVerticesIndx.size(); iVert++)
		{
			bigSegmentColors.row(regionVerticesIndx[iVert]) << colorIt->second.coeff(0), colorIt->second.coeff(1), colorIt->second.coeff(2), 1;
		}
		colorIt++;
	}


	std::set<int> ixCurrSegBigNeigh;
	std::set<int> ixCurrSegSmallNeigh;
	std::map<int, int> segSrc2segDst; // The src seg that will be merge to dst seg

	std::cout << "Start the merge process" << std::endl;

	while (!smallSegments.empty())
	{
		std::cout << smallSegments.size() << " small segments remained to be merged" << std::endl;
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
				if (bigSegments.count(iNeighboorSeg) > 0)
				{
					ixCurrSegBigNeigh.insert(iNeighboorSeg);
				}
				else
				{
					ixCurrSegSmallNeigh.insert(iNeighboorSeg);
				}
			}

			/*
				Find the most similar in the big segments nieghbors
			*/

			
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


			/*for (std::set<int>::iterator itSmallNeigh = ixCurrSegSmallNeigh.begin(); itSmallNeigh != ixCurrSegSmallNeigh.end(); itSmallNeigh++)
			{

				double currDiff = abs(segmentsAvgCurvedness.at(*itSmallNeigh) - segmentsAvgCurvedness.at(smallSegIt->first));
				if (currDiff < minCurvDiff)
				{
					iBigSegMostSim = *itSmallNeigh;
					minCurvDiff = currDiff;
				}
			}*/

			int iBigSegMostSim = -1;
			double maxInnerProd = -2; // The vectors should be normalized
			//double minCurvDiff = 99999;

			for (std::set<int>::iterator itBigNeigh = ixCurrSegBigNeigh.begin(); itBigNeigh != ixCurrSegBigNeigh.end(); itBigNeigh++)
			{

				double currInnerProd = segmentsAvgNormal.at(*itBigNeigh).dot(segmentsAvgNormal.at(smallSegIt->first));
				if (currInnerProd > maxInnerProd)
				{
					iBigSegMostSim = *itBigNeigh;
					maxInnerProd = currInnerProd;
				}
			}

			//for (std::set<int>::iterator itSmallNeigh = ixCurrSegSmallNeigh.begin(); itSmallNeigh != ixCurrSegSmallNeigh.end(); itSmallNeigh++)
			//{

			//	double currInnerProd = segmentsAvgNormal.at(*itSmallNeigh).dot(segmentsAvgNormal.at(smallSegIt->first));
			//	if (currInnerProd > maxInnerProd)
			//	{
			//		iBigSegMostSim = *itSmallNeigh;
			//		maxInnerProd = currInnerProd;
			//	}
			//}

			if (iBigSegMostSim != -1)
			{
				segSrc2segDst.insert({ smallSegIt->first,iBigSegMostSim });
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

			for (int iVert = 0; iVert < smallSegments.at(iSrcSeg).piece_vertices_index_.size(); iVert++)
			{
				vertIndex2SegIndex[iVert] = iDstSeg;
			}

			if (bigSegments.count(iDstSeg) > 0)
			{
				bigSegments.at(iDstSeg).piece_vertices_index_.insert(
					bigSegments.at(iDstSeg).piece_vertices_index_.end(),
					smallSegments.at(iSrcSeg).piece_vertices_index_.begin(),
					smallSegments.at(iSrcSeg).piece_vertices_index_.end()
				);
			}
			else {
				/*smallSegments.at(iDstSeg).piece_vertices_index_.insert(
					smallSegments.at(iDstSeg).piece_vertices_index_.end(),
					smallSegments.at(iSrcSeg).piece_vertices_index_.begin(),
					smallSegments.at(iSrcSeg).piece_vertices_index_.end()
				);*/
			}
			

			oRegionOutsideBoundaryVerticesList[iDstSeg].insert(
				oRegionOutsideBoundaryVerticesList[iDstSeg].end(),
				oRegionOutsideBoundaryVerticesList[iSrcSeg].begin(),
				oRegionOutsideBoundaryVerticesList[iSrcSeg].end()
			);

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