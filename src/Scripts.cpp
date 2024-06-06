#include "Scripts.h"
#include "Visualization.h"
#include "ScriptsUtils.h"
#include "SpectralClustering.h"


Segment segment_intact_surface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave, bool isVisualizer)
{

	std::cout << "**************** segment_intact_surface **************** " << std::endl;

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	std::map<int, Segment*> smallSegments;
	std::map<int, Segment*> bigSegments;
	double fracture = intactSimilarityFracture; // 0.65; //0.65;
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
		
		/*if (bigSegments.size() == 1)
		{
			fracture = fracture - 0.12;
			continue;
		}*/
		

		intactIndex = 0;//-1;
		double minMeanCurvedness = 9999999;

		for(auto seg = bigSegments.begin(); seg!= bigSegments.end(); ++seg)
		{
			double segCurvedness = 0;
			//auto piece_vertices_index_ = seg.piece_vertices_index_;
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
				intactIndex = seg->first; //k;
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

		/*if (isSaveAll)
		{
			int ii = 0;

			for (auto seg = bigSegments.begin(); seg != bigSegments.end(); ++seg) 
			{
				seg->second->saveAsObj(fragment.m_FolderPath + "\\" + fragment.m_Name + "_intact.obj");

			}
		}*/

		segments[intactIndex].saveAsObj(fragment.m_FolderPath + "\\" + fragment.m_Name+ "_intact.obj" );
		std::cout << "Write successfully the output to path " << fragment.m_FolderPath << std::endl;
	}

	//bool isVisualizer = true; // Make this function param when refactor
	if (isVisualizer)
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

void segment_opposite_surface(ObjFragment& fragment, double intactNormalStdThershold,  double intactSimilarityFracture, bool isSave, bool isVisualizer)
{
	Segment intactSegment = segment_intact_surface(fragment, intactNormalStdThershold,intactSimilarityFracture,isSave,false);

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
	double ABUTTING_PERCENTAGE = 0.1; //0.05; // we have a bug on updating the abutting so the thresh is low

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

		/// This if is just beacue merge is expecting map....
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

	if (isVisualizer)
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



void segment_sidewalls_surface(ObjFragment& fragment,double intactNormalStdThershold, double intactSimilarityFracture, bool isSave, bool isVisualizer)
{
	Segment intactSegment = segment_intact_surface(fragment, intactNormalStdThershold,intactSimilarityFracture,isSave,false);

	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	std::map<int, Segment*> smallSegments;
	std::map<int, Segment*> bigSegments;
	double fracture = 0.5;//0.5;
	double minBigSegPercSize = (double)1 / 8000; //0.00025;
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
	
	if (isVisualizer)
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

void colorSmooth(ObjFragment& fragment, bool isSave, bool isVisualizer)
{
	/*

		fragment is essentially could be segment ....
	*/
	double fracture = 0.2; //0.3; //0.65;
	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	std::map<int, int> vertIndex2SegIndex;
	std::vector<Segment> segments;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	double simThresh = fragment.getSimilarThreshByPos(fracture);
	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);
	initSegments(segments, vertIndex2SegIndex, oRegionsList, oRegionOutsideBoundaryVerticesList, fragmentSize, fragment);

	std::map<int, Eigen::RowVector3d> meshColors;
	generateRandomColors(meshColors, oRegionsList.size());
	auto colorIt = meshColors.begin();
	/*std::map<int, Eigen::RowVector3d> tfiraBlackWhite;
	tfiraBlackWhite.insert({ 0, { 0,0,0 } });
	tfiraBlackWhite.insert({ 1, { 0.5,0.5,0.5 } });*/
	auto smoothColorIt = meshColors.begin();
	auto unSmoothColorIt = std::next(meshColors.begin(),45);

	double smooth_threshold = fragment.getSimilarThreshByPos(0.25);
	Eigen::MatrixXd segment2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::MatrixXd smooth2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::MatrixXd frag2UnionColor = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	int iSeg = 0;

	for (Segment& seg : segments)
	{
		/*
		
		NOTE THAT THIS IS HARDCODED, YOU HAVE THE SAME FUNCTIONALITY IN SEGMENT INTACT
		*/

		double segCurvedness = 0;
		auto& piece_vertices_index_ = seg.piece_vertices_index_;

		for (int verIndex : piece_vertices_index_)
		{
			segCurvedness += fragment.m_NormedMeshCurvedness[verIndex];
		}

		double segAvgCur = segCurvedness / piece_vertices_index_.size();

		if (segAvgCur < smooth_threshold)
		{
			colorFragSingleSeg(smooth2Colors, seg, smoothColorIt);
		}
		else
		{
			colorFragSingleSeg(smooth2Colors, seg, unSmoothColorIt);
		}

		colorFragSingleSeg(segment2Colors, seg, colorIt++);
		colorFragSingleSeg(frag2UnionColor, seg, smoothColorIt);
	}

	Eigen::MatrixXd vertsCurvdness2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	for (int i = 0; i < fragment.m_NormedMeshCurvedness.size(); i++)
	{
		double val = fragment.m_NormedMeshCurvedness(i);
		vertsCurvdness2Colors.row(i) << val, 0, 0, 1;
	}


	Eigen::MatrixXd vertsCoordsNormalized = fragment.m_Vertices;

	// To enlarge the variance
	// var(ax) = a**2var(x)
	double factor = 99999999;
	for (int i = 0; i < vertsCoordsNormalized.rows(); i++)
	{

		vertsCoordsNormalized.row(i) << vertsCoordsNormalized.coeff(i, 0) * factor, vertsCoordsNormalized.coeff(i, 1) * factor, vertsCoordsNormalized.coeff(i, 2) * factor;
	}

	for (int i = 0; i < 3; i++)
	{
		vertsCoordsNormalized.col(i).normalize();
	}

	Eigen::MatrixXd fragVerts2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::MatrixXd vertsHieghtworldAxis2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	for (int i = 0; i < fragment.m_NormedMeshCurvedness.size(); i++)
	{
		double val = fragment.m_NormedMeshCurvedness(i);// * 2;
		double heightBias = vertsCoordsNormalized.coeff(i,2) * 100; // double 10000 so the numbers will be 0.x and not 0.000x. we have here an assumption that the opposite surface is rotated..
		//double heightBias = fragment.m_Vertices.coeff(i, 2) * 100; // double 10000 so the numbers will be 0.x and not 0.000x. we have here an assumption that the opposite surface is rotated..
		//fragVerts2Colors.row(i) << val, val, val,val;
		//fragVerts2Colors.row(i) << val + vertsCoordsNormalized.coeff(i,2), 0, 0, val;
		//fragVerts2Colors.row(i) << val + vertsCoordsNormalized.coeff(i,2)*10000, 0, 0, val;
		//heightBias = heightBias * heightBias;
		//fragVerts2Colors.row(i) << heightBias*10, 0, 0, val*heightBias;
		//fragVerts2Colors.row(i) << heightBias, 0, 0, val+heightBias* heightBias* heightBias;
		//fragVerts2Colors.row(i) << heightBias, 0, 0, val - heightBias;//* heightBias* heightBias;


		fragVerts2Colors.row(i) << val - heightBias, 0, 0, 1;//* heightBias* heightBias;
		//fragVerts2Colors.row(i) << val - std::exp(-(1-heightBias + 0.5)), 0, 0, 1;//* heightBias* heightBias;
		vertsHieghtworldAxis2Colors.row(i) << heightBias, 0, 0, 1;
	}






	//https://forum.kde.org/viewtopic.php?f=9&t=110265
	Eigen::MatrixXd fragAligned = fragment.m_Vertices.rowwise() - fragment.m_Vertices.colwise().mean();

	// we can directly take SVD
	Eigen::JacobiSVD<Eigen::MatrixXd> svd(fragAligned, Eigen::ComputeThinV);

	Eigen::MatrixXd W = svd.matrixV().leftCols(3);

	// then to project the data:
	Eigen::MatrixXd projected = fragAligned * W; 

	Eigen::MatrixXd fragVertsInSvdcol02Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::VectorXd height_0 = projected.col(0); // hardcoded for RPF_438 group 52 , this is the height of the manifold
	//heightNormed = heightNormed.normalized() * 100;
	fragVertsInSvdcol02Colors.col(0) = height_0;
	fragVertsInSvdcol02Colors.col(3) = Eigen::VectorXd::Ones(height_0.rows());
	
	Eigen::MatrixXd fragVertsInSvdcol12Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::VectorXd height_1 = projected.col(1);
	//heightNormed  = heightNormed.normalized() * -100;
	fragVertsInSvdcol12Colors.col(0) = height_1;
	fragVertsInSvdcol12Colors.col(3) = Eigen::VectorXd::Ones(height_1.rows());

	Eigen::MatrixXd fragVertsInSvdcol22Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::VectorXd height_2 = projected.col(2);
	//heightNormed  = heightNormed.normalized() * -100;
	fragVertsInSvdcol22Colors.col(0) = height_2;
	fragVertsInSvdcol22Colors.col(3) = Eigen::VectorXd::Ones(height_2.rows());
	
	Eigen::MatrixXd fragVertsInSvdcol2Normed2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);
	Eigen::VectorXd heightNormed_2 = projected.col(2);
	heightNormed_2 = heightNormed_2.normalized() * 10;//0.xx and not 0.0xx
	fragVertsInSvdcol2Normed2Colors.col(0) = heightNormed_2;
	fragVertsInSvdcol2Normed2Colors.col(3) = Eigen::VectorXd::Ones(heightNormed_2.rows());
	
	for (int i = 0; i < heightNormed_2.size(); i++)
	{
		if (sigmiod(heightNormed_2.coeff(i))>0.5)
		{
			fragVertsInSvdcol2Normed2Colors(i, 0) = 1;
		}
		else
		{
			fragVertsInSvdcol2Normed2Colors(i, 0) = 0;
		}
	}


	// for RPF_438 group 52
	Eigen::MatrixXd fragVertsEigenValueHieght2Colors = Eigen::MatrixXd::Zero(fragment.m_Vertices.rows(), 4);


	std::cout << "max normed curv" << fragment.m_NormedMeshCurvedness.maxCoeff() << std::endl;

	for (int i = 0; i < fragment.m_NormedMeshCurvedness.size(); i++)
	{
		double curv = fragment.m_NormedMeshCurvedness(i);// * 2;
		curv *= 1.2;
		double height = heightNormed_2(i); /// Assume it has the lowest variance
		double downgrade = std::exp(-(height+0.42)); //0.08 * 1/(height);
		//double downgrade = std::exp(-(height*1.0001)); //0.08 * 1/(height);
		double r = curv - downgrade;
		
		if (sigmiod(r) >= 0.525)
		{
			r = 1;
		}
		else
		{
			r = 0;
		}

		fragVertsEigenValueHieght2Colors.row(i) << r, 0, 0,1 ;//* heightBias* heightBias;
	}

	Eigen::MatrixXd fragNormalAligned = fragment.m_Normals.rowwise() - fragment.m_Normals.colwise().mean();

	// we can directly take SVD
	Eigen::JacobiSVD<Eigen::MatrixXd> svdNormal(fragNormalAligned, Eigen::ComputeThinV);

	Eigen::MatrixXd Wnormal = svdNormal.matrixV().leftCols(3);

	// then to project the data:
	Eigen::MatrixXd normalProjected = fragNormalAligned * Wnormal;


	Eigen::MatrixXd normalSvdcol02Colors = Eigen::MatrixXd::Zero(fragment.m_Normals.rows(), 4);
	Eigen::VectorXd normal_0 = normalProjected.col(0); // hardcoded for RPF_438 group 52 , this is the height of the manifold
	//heightNormed = heightNormed.normalized() * 100;
	normalSvdcol02Colors.col(0) = normal_0;
	normalSvdcol02Colors.col(3) = Eigen::VectorXd::Ones(normal_0.rows());

	Eigen::MatrixXd normalSvdcol12Colors = Eigen::MatrixXd::Zero(fragment.m_Normals.rows(), 4);
	Eigen::VectorXd normal_1 = normalProjected.col(1);
	//heightNormed  = heightNormed.normalized() * -100;
	normalSvdcol12Colors.col(0) = normal_1;
	normalSvdcol12Colors.col(3) = Eigen::VectorXd::Ones(normal_1.rows());

	Eigen::MatrixXd normalSvdcol22Colors = Eigen::MatrixXd::Zero(fragment.m_Normals.rows(), 4);
	Eigen::VectorXd normal_2 = normalProjected.col(2);
	//heightNormed  = heightNormed.normalized() * -100;
	normalSvdcol22Colors.col(0) = normal_2;
	normalSvdcol22Colors.col(3) = Eigen::VectorXd::Ones(normal_2.rows());

	if (isVisualizer)
	{
		Visualizer visualizer;

		visualizer.appendMesh(fragment.m_Vertices, fragment.m_Faces, meshColors[0]);
		visualizer.m_Viewer.callback_key_down =
			[&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) ->bool
		{

			switch (key) {
			case '1':
				visualizer.m_Viewer.data().set_colors(segment2Colors);//meshColors[0]
				std::cout << "Pressed 1" << std::endl;
				break;

			case '2':
				visualizer.m_Viewer.data().set_colors(smooth2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

			case '3':
				visualizer.m_Viewer.data().set_colors(frag2UnionColor);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
				return false;
			
			case '4':
				visualizer.m_Viewer.data().set_colors(vertsCurvdness2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

			case '5':
				visualizer.m_Viewer.data().set_colors(fragVerts2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			
			case '6':
				visualizer.m_Viewer.data().set_colors(vertsHieghtworldAxis2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

			case '7':
				visualizer.m_Viewer.data().set_colors(fragVertsInSvdcol02Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			case '8':
				visualizer.m_Viewer.data().set_colors(fragVertsInSvdcol12Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

			case '9':
				visualizer.m_Viewer.data().set_colors(fragVertsInSvdcol22Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			case '0':
				visualizer.m_Viewer.data().set_colors(fragVertsInSvdcol2Normed2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			case '-':
				visualizer.m_Viewer.data().set_colors(fragVertsEigenValueHieght2Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			case 'Q':
				visualizer.m_Viewer.data().set_colors(normalSvdcol02Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			case 'W':
				visualizer.m_Viewer.data().set_colors(normalSvdcol12Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;
			case 'E':
				visualizer.m_Viewer.data().set_colors(normalSvdcol22Colors);//meshColors[0]
				std::cout << "Pressed 2" << std::endl;
				break;

				return false;
			};

			
	
		};
		visualizer.launch();
	}


	if (isSave)
	{
		std::string mtlName = fragment.m_Name + "_bambooMap.mtl";
		std::string mtlPath = fragment.m_FolderPath + "\\" + mtlName;
		std::string imgPath = "platte.jpg";
		std::string materialName = "material_0";
		saveMtlFile(mtlPath, imgPath);

		
		int knownColors = 2;
		Eigen::MatrixXd textureCoords(knownColors*3,2);
		
		// red color
		textureCoords.row(0) << 1.000000, 0.333333;
		textureCoords.row(1) << 1.000001, 0.333335;
		textureCoords.row(2) << 1.000002, 0.333333;
		
		// white color
		textureCoords.row(3) << 1.000000, 1.000000;
		textureCoords.row(4) << 1.000000, 1.000005;
		textureCoords.row(5) << 0.999999, 1.000003;

		//double heightMean = fragVertsInSvdcol22Colors.col(0).mean();

		std::vector<double> heights_;
		for (int64_t i = 0; i < fragVertsInSvdcol22Colors.col(0).rows(); ++i)
			heights_.push_back(fragVertsInSvdcol22Colors(i,0));

		std::sort(heights_.begin(), heights_.end());
		const auto median_it = heights_.begin() + heights_.size() * 0.70;//0.65 good for fine segmentation  //0.5 good for detecting intact vs fractured when removing region merging 
		std::nth_element(heights_.begin(), median_it, heights_.end());
		double heightMean = *median_it;
		
		Eigen::MatrixXi face2TextureCoords(fragment.m_Faces.rows(),3);

		for (int i = 0; i < fragment.m_Faces.rows(); i+=1)
		{
			int vert1 = fragment.m_Faces(i, 0);
			int vert2 = fragment.m_Faces(i, 1);
			int vert3 = fragment.m_Faces(i, 2);
			double h = (fragVertsInSvdcol22Colors(vert1, 0) + fragVertsInSvdcol22Colors(vert2, 0) + fragVertsInSvdcol22Colors(vert3, 0)) / 3;
			
			if (h > heightMean)
			{
				face2TextureCoords.row(i) << 0, 1, 2;
			}else
			{
				face2TextureCoords.row(i) << 3, 4, 5;
			}
		}
		
		std::string objPath = fragment.m_FolderPath + "\\" + fragment.m_Name + "_bambooMap.obj";
		saveObjFile(objPath, mtlName,
			fragment.m_Vertices, fragment.m_Faces, fragment.m_Normals,
			fragment.m_Faces2Normals, textureCoords, face2TextureCoords);

	}
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