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

	double minSegPercSize = 0.005; // This need to get as a parameter to the script
	std::map<int, Segment> smallSegments;
	std::map<int, Segment> bigSegments;
	std::map<int, double> segmentsAvgCurvedness; // We avg? maybe std is a parameter we should consider
	std::map<int, int> vertIndex2SegIndex;
	double fragmentSize = static_cast<double>(fragment.m_Vertices.rows());

	for (int iSeg = 0; iSeg < oRegionsList.size(); iSeg++)
	{
		auto segVertsIndexes = oRegionsList[iSeg];
		Segment currSeg(segVertsIndexes);
		segmentsAvgCurvedness.insert({ iSeg,fragment.localCurvedness(segVertsIndexes) });
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