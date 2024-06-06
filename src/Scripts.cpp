#include "Scripts.h"
#include "Visualization.h"
#include "ScriptsUtils.h"


Segment segment_intact_surface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture,
	bool isSave, bool isDebug)
{
	std::cout << "**************** segment_intact_surface **************** " << std::endl;

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	std::map<int, Segment*> smallSegments;
	std::map<int, Segment*> bigSegments;
	double fracture = intactSimilarityFracture; 
	double minBigSegPercSize = 0.05;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());
	double MAX_NUM_TRIALS = 6;

	bool isSegmented = false;
	int intactIndex = -1;
	int nTrials = 0;

	while (!isSegmented)
	{
		if (nTrials >= MAX_NUM_TRIALS)
		{
			std::cout << "The number of trials exceeded the number of maximum allocated (" << MAX_NUM_TRIALS << ")" << std::endl;
			std::cout << "Please rerun mannually..exiting" << std::endl;
			std::exit(1);
		}

		if (fracture > 1 || fracture <=0)
		{
			std::cout << "fracrure value is " << fracture << " and it should be between 0 to 1" << std::endl;
			std::cout << "Please rerun " << fragment.m_Name << " mannually..exiting" << std::endl;
			std::exit(1);
		}

		oRegionsList.clear();
		oRegionOutsideBoundaryVerticesList.clear();
		vertIndex2SegIndex.clear();
		segments.clear();
		smallSegments.clear();
		bigSegments.clear();
		nTrials += 1;

		double simThresh = fragment.getSimilarThreshByPos(fracture);
		std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
		fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);

		initSegments(segments, vertIndex2SegIndex, oRegionsList, oRegionOutsideBoundaryVerticesList, fragmentSize, fragment);
		sortToSmallAndBigSegments(smallSegments, bigSegments, segments, minBigSegPercSize);

		std::cout << "In trial: " << nTrials << " Found " << bigSegments.size() << "big Segments" << std::endl;

		if (bigSegments.size() == 0)
		{
			fracture = fracture + 0.05;
			continue;
		}
		
		intactIndex = 0;
		double minMeanCurvedness = 9999999;

		for(auto seg = bigSegments.begin(); seg!= bigSegments.end(); ++seg)
		{
			double segCurvedness = 0;
			auto& piece_vertices_index_ = seg->second->piece_vertices_index_;

			for (int verIndex : piece_vertices_index_)
			{
				segCurvedness += fragment.m_NormedMeshCurvedness[verIndex];
			}
			
			seg->second->loadNormedNormals();
			Eigen::Vector3d avgNormal_debug = calcAvg(seg->second->m_NormedNormals);
			auto y = avgNormal_debug.coeff(1);

			// checking the intact turns its face up
			if (y < 0)
			{
				continue;
			}

			double segAvgCur = segCurvedness / piece_vertices_index_.size();

			if (segAvgCur < minMeanCurvedness)
			{
				minMeanCurvedness = segAvgCur;
				intactIndex = seg->first; 
			}
		}

		Segment& intactSurface = segments[intactIndex];
		//intactSurface.loadNormedNormals();
		Eigen::Vector3d avgNormal = calcAvg(intactSurface.m_NormedNormals);
		Eigen::Vector3d stdNormal = calcVariance(intactSurface.m_NormedNormals, avgNormal).array().sqrt();
		double l2 = stdNormal.norm();

		if (l2 < intactNormalStdThershold)
		{
			isSegmented = true;
		}
		else {
			std::cout << "the std of the normals is " << l2 <<" , it seems to be to much big - continue to segment" << std::endl;
			fracture = fracture - 0.08;
		}
	}

	segments[intactIndex].loadBasicData();
	
	if (isSave)
	{
		segments[intactIndex].saveAsObj(fragment.m_FolderPath + "\\" + fragment.m_Name+ "_intact.obj" );
		std::cout << "Write successfully the output to path " << fragment.m_FolderPath << std::endl;
	}

	if (isDebug)
	{
		std::map<int, Eigen::RowVector3d> meshColors;
		generateRandomColors(meshColors, segments.size());
		
		auto colorIt = meshColors.begin();
		Eigen::MatrixXd segment2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
		colorFrag(segment2Colors, smallSegments, colorIt);
		colorFrag(segment2Colors, bigSegments, colorIt);

		colorIt = meshColors.begin();
		Eigen::MatrixXd intact2Colors = Eigen::MatrixXd::Ones(fragment.m_Vertices.rows(), 4);
		std::map<int, Segment*> tmpContainer;
		tmpContainer.insert({intactIndex,&segments[intactIndex]});
		colorFrag(intact2Colors, tmpContainer, colorIt);

		colorIt = meshColors.begin();
		Eigen::MatrixXd bigSegments2Colors = Eigen::MatrixXd::Ones(fragment.m_Vertices.rows(), 4);
		colorFrag(bigSegments2Colors, bigSegments, colorIt);

		Visualizer visualizer;
		visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, fragment.m_Colors);
		visualizer.m_Viewer.callback_key_down =
			[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) ->bool
		{

			switch (key) {
			case '1':
				visualizer.m_Viewer.data().set_colors(segment2Colors);//meshColors[0]
				std::cout << "Pressed 1" << std::endl;
				break;
			case '2':
				visualizer.m_Viewer.data().set_colors(intact2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

			case '3':
				visualizer.m_Viewer.data().set_colors(bigSegments2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

				return false;
			};

		};
		visualizer.launch();
	}

	std::cout << "Return the intact surface" << std::endl;
	return segments[intactIndex];
}

void segment_opposite_surface(ObjFragment& fragment, double intactNormalStdThershold,  double intactSimilarityFracture, 
	bool isSave, bool isDebug)
{
	Segment intactSegment = segment_intact_surface(fragment, intactNormalStdThershold,intactSimilarityFracture,false,false);

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	std::map<int, Segment*> smallSegments;
	std::map<int, Segment*> bigSegments;
	double fracture = 0.5;
	double minBigSegPercSize = 0.0025;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	double simThresh = fragment.getSimilarThreshByPos(fracture);
	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);
	initSegments(segments, vertIndex2SegIndex, oRegionsList, oRegionOutsideBoundaryVerticesList, fragmentSize, fragment);
	sortToSmallAndBigSegments(smallSegments, bigSegments, segments, minBigSegPercSize);

	std::map<int, Eigen::RowVector3d> meshColors;
	generateRandomColors(meshColors, oRegionsList.size());
	

	Eigen::MatrixXd rawBigSegmentColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(rawBigSegmentColors, bigSegments, meshColors.begin());
	Eigen::MatrixXd rawSegmentsColors = rawBigSegmentColors;
	colorFrag(rawSegmentsColors, smallSegments, std::next(meshColors.begin(), bigSegments.size()));
	//mergeSmall2BigSegments(smallSegments, bigSegments, vertIndex2SegIndex);
	Eigen::MatrixXd afterBig2SmallColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	//colorFrag(afterBig2SmallColors, bigSegments, meshColors.begin());

	Eigen::Vector3d intactAvgNormal = calcAvg(intactSegment.m_NormedNormals);
	std::map<int, Segment*> oppositeSegSeeds;
	double EPSILONERR = 0.2;

	for (int i=0;i< segments.size();++i)
	{
		segments[i].loadNormedNormals();
		Eigen::Vector3d currSegAvg = calcAvg(segments[i].m_NormedNormals);

		if (currSegAvg.dot(intactAvgNormal) < -(1 - EPSILONERR))
		{
			oppositeSegSeeds.insert({ i,&segments[i]});
		}
	}

	Eigen::MatrixXd oppSegmentSeedColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(oppSegmentSeedColors, oppositeSegSeeds, meshColors.begin());
	auto oppSegIt = oppositeSegSeeds.begin();
	std::vector<int> segMergedIndexes;

	for (std::map<int, Segment*>::iterator it = std::next(oppositeSegSeeds.begin(),1); it != oppositeSegSeeds.end(); it++)
	{
		merge(it->first, oppSegIt->first, vertIndex2SegIndex, oppositeSegSeeds, oppositeSegSeeds);
		segMergedIndexes.push_back(it->first);
	}

	for (int i=0; i < segMergedIndexes.size(); ++i)
	{

		oppositeSegSeeds.erase(segMergedIndexes[i]);
		if (smallSegments.count(segMergedIndexes[i]) > 0)
		{
			smallSegments.erase(segMergedIndexes[i]);
		}
		else
		{
			bigSegments.erase(segMergedIndexes[i]);
		}
	}
	
	auto colorIt = meshColors.begin();
	Eigen::MatrixXd oppSegment2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(oppSegment2Colors, oppositeSegSeeds, colorIt);

	std::map<int, int> iNeigh2Count;
	oppSegIt->second->findNeighbors(iNeigh2Count, vertIndex2SegIndex, oppSegIt->first);
	segMergedIndexes.clear();
	double ABUTTING_PERCENTAGE = 0.1; //0.05; 

	for (auto it = iNeigh2Count.begin(); it != iNeigh2Count.end(); it++)
	{
		int segBoundarySize = segments[it->first].m_OutsideBoundaryVertsIndexes.size(); 
		double abuttingPercentage = it->second / static_cast<double>(segBoundarySize);
		if (abuttingPercentage > ABUTTING_PERCENTAGE)
		{
			segMergedIndexes.push_back(it->first);

			// For debugging 
			colorFragSingleSeg(oppSegment2Colors, segments[it->first], colorIt);
			++colorIt;

		}

	}

	

	for (int i = 0; i < segMergedIndexes.size(); ++i)
	{

		if (smallSegments.count(segMergedIndexes[i]) > 0)
		{
			merge(segMergedIndexes[i], oppSegIt->first, vertIndex2SegIndex, smallSegments, oppositeSegSeeds);
			smallSegments.erase(segMergedIndexes[i]);
		}
		else
		{
			merge(segMergedIndexes[i], oppSegIt->first, vertIndex2SegIndex, bigSegments, oppositeSegSeeds);
			bigSegments.erase(segMergedIndexes[i]);
		}
	}

	if (isSave)
	{
	
		if (oppositeSegSeeds.size() !=1)
		{
			std::cout << "Expected oppositeSegSeeds to contain one element but it is empty" << std::endl;
		}
		else
		{
			auto finalSeg = oppositeSegSeeds.begin()->second;
			finalSeg->loadBasicData();
			finalSeg->saveAsObj(fragment.m_FolderPath + "\\" + fragment.m_Name + "_opposite.obj");
			std::cout << "Write opposite surface successfully the output to path " << fragment.m_FolderPath << std::endl;
		}
	}

	if (isDebug)
	{


		Visualizer visualizer;
		visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, meshColors[0]);
		visualizer.m_Viewer.callback_key_down =
			[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) ->bool
		{

			switch (key) {
			case '1':
				visualizer.m_Viewer.data().set_colors(fragment.m_Colors);//meshColors[0]
				std::cout << "Pressed 1" << std::endl;
				break;

			case '2':
				visualizer.m_Viewer.data().set_colors(rawSegmentsColors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

			case '3':
				visualizer.m_Viewer.data().set_colors(rawBigSegmentColors);
				std::cout << "Pressed 3, present only the big segments" << std::endl;
				break;
			case '4':
				visualizer.m_Viewer.data().set_colors(afterBig2SmallColors);
				std::cout << "Pressed 4, present only the big segments" << std::endl;
				break;
			case '5':
				visualizer.m_Viewer.data().set_colors(oppSegmentSeedColors);
				std::cout << "Pressed 5, present only the big segments" << std::endl;
				break;

			case '6':
				visualizer.m_Viewer.data().set_colors(oppSegment2Colors);
				std::cout << "Pressed 6, present only the big segments" << std::endl;
				break;

				return false;
			};

		};
		visualizer.launch();
	}
}



void segment_sidewalls_surface(ObjFragment& fragment,double intactNormalStdThershold, double intactSimilarityFracture,
	bool isSave, bool isDebug)
{
	Segment intactSegment = segment_intact_surface(fragment, intactNormalStdThershold,intactSimilarityFracture,false,false);

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	std::map<int, Segment*> smallSegments;
	std::map<int, Segment*> bigSegments;
	double fracture = 0.5;
	double minBigSegPercSize = (double)1 / 8000; 
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	double simThresh = fragment.getSimilarThreshByPos(fracture);
	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);
	initSegments(segments, vertIndex2SegIndex, oRegionsList, oRegionOutsideBoundaryVerticesList, fragmentSize, fragment);
	sortToSmallAndBigSegments(smallSegments, bigSegments, segments, minBigSegPercSize);
	std::map<int, Eigen::RowVector3d> meshColors;
	generateRandomColors(meshColors, oRegionsList.size());
	

	Eigen::MatrixXd rawBigSegmentColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(rawBigSegmentColors, bigSegments, meshColors.begin());
	Eigen::MatrixXd rawSegmentsColors = rawBigSegmentColors;
	colorFrag(rawSegmentsColors, smallSegments, std::next(meshColors.begin(), bigSegments.size()));
	mergeSmall2BigSegments(smallSegments, bigSegments, vertIndex2SegIndex);
	Eigen::MatrixXd afterBig2SmallColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(afterBig2SmallColors, bigSegments, meshColors.begin());

	Eigen::Vector3d intactAvgNormal = calcAvg(intactSegment.m_NormedNormals);
	std::map<int, Segment*> sideWallSegSeeds;
	double EPSILON_ERR = 0.5;

	for (auto &it : bigSegments)
	{
		auto& seg = it.second;
		auto& ix = it.first;
		seg->loadNormedNormals();
		seg->avgNormedNormal = calcAvg(seg->m_NormedNormals);

		double dotProduct = seg->avgNormedNormal.dot(intactAvgNormal);

		if ( abs(dotProduct) < EPSILON_ERR)
		{
			sideWallSegSeeds.insert(it);
		}
	}

	Eigen::MatrixXd wallSegmentSeedColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(wallSegmentSeedColors, sideWallSegSeeds, meshColors.begin());

	
	/*double EPSILON_ERR2 = 0.2;
	std::map<int, int> segSrc2segDst;
	std::vector<int> potentialSeeds;
	std::set<int> mergedAlreadySegIndexes;
	int nLastSize = -1;

	while(mergedAlreadySegIndexes.size()!=nLastSize)
	{
		nLastSize = mergedAlreadySegIndexes.size();

		for (auto& itSegSeed : sideWallSegSeeds)
		{
			std::map<int, int> iNeigh2Count;
			auto& currSeg = itSegSeed.second;
			int currSegSeedIndex = itSegSeed.first;
			currSeg->findNeighbors(iNeigh2Count, vertIndex2SegIndex, itSegSeed.first);


			for (auto& nbr : iNeigh2Count)
			{
				if (sideWallSegSeeds.find(nbr.first) != sideWallSegSeeds.end())
				{
					if (mergedAlreadySegIndexes.find(nbr.first) == mergedAlreadySegIndexes.end())
					{
						auto& nbrSeg = sideWallSegSeeds.at(nbr.first);
						double dotProduct = currSeg->avgNormedNormal.dot(nbrSeg->avgNormedNormal);

						if (1 - abs(dotProduct) < EPSILON_ERR2)
						{
							segSrc2segDst.insert({ currSegSeedIndex ,nbr.first });
						}

					}
				}
			}

			for (std::map<int, int>::iterator itSegSrc2Dst = segSrc2segDst.begin(); itSegSrc2Dst != segSrc2segDst.end(); ++itSegSrc2Dst)
			{
				int iSrcSeg = itSegSrc2Dst->first;
				int iDstSeg = itSegSrc2Dst->second;
				merge(iSrcSeg, iDstSeg, vertIndex2SegIndex, sideWallSegSeeds, sideWallSegSeeds);
				currSeg->loadNormedNormals();
				currSeg->avgNormedNormal = calcAvg(currSeg->m_NormedNormals);
			}

			for (std::map<int, int>::iterator itSrc2Dst = segSrc2segDst.begin(); itSrc2Dst != segSrc2segDst.end(); ++itSrc2Dst)
			{
				mergedAlreadySegIndexes.emplace(itSrc2Dst->first);
			}

			segSrc2segDst.clear();

		}
		
	} 

	*/

	
	auto firstSeed = sideWallSegSeeds.begin();

	for(auto itSegSeed = std::next(sideWallSegSeeds.begin(),1); itSegSeed!=sideWallSegSeeds.end() ;)
	{
		merge(itSegSeed->second, firstSeed->second, firstSeed->first,vertIndex2SegIndex);
		itSegSeed = sideWallSegSeeds.erase(itSegSeed);
	}

	if (isSave)
	{
		if (sideWallSegSeeds.size() != 1)
		{
			std::cout << "Expected sideWallSegSeeds to contain one element but it is empty" << std::endl;
		}
		else
		{
			auto finalSeg = sideWallSegSeeds.begin()->second;
			finalSeg->loadBasicData();
			finalSeg->saveAsObj(fragment.m_FolderPath + "\\" + fragment.m_Name + "_sideWalls.obj");
			std::cout << "Write side walls surface successfully the output to path " << fragment.m_FolderPath << std::endl;
		}
	}

	Eigen::MatrixXd wallSegmentMergeColors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	colorFrag(wallSegmentMergeColors, sideWallSegSeeds, meshColors.begin());
	
	if (isDebug)
	{
		Visualizer visualizer;

		visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, meshColors[0]);
		visualizer.m_Viewer.callback_key_down =
			[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) ->bool
		{

			switch (key) {
			case '1':
				visualizer.m_Viewer.data().set_colors(fragment.m_Colors);
				std::cout << "Pressed 1" << std::endl;
				break;

			case '2':
				visualizer.m_Viewer.data().set_colors(rawSegmentsColors);
				std::cout << "Pressed 2" << std::endl;
				break;

			case '3':
				visualizer.m_Viewer.data().set_colors(rawBigSegmentColors);
				std::cout << "Pressed 3, present only the big segments" << std::endl;
				break;
			case '4':
				visualizer.m_Viewer.data().set_colors(afterBig2SmallColors);
				std::cout << "Pressed 4, present only the big segments" << std::endl;
				break;
			case '5':
				visualizer.m_Viewer.data().set_colors(wallSegmentSeedColors);
				std::cout << "Pressed 5, present only the big segments" << std::endl;
				break;
			case '6':
				visualizer.m_Viewer.data().set_colors(wallSegmentMergeColors);
				std::cout << "Pressed 5, present only the big segments" << std::endl;
				break;

				return false;
			};

		};
		visualizer.launch();
	}
}





