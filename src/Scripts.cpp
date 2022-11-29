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
	double fracture = 0.8;
	double simThresh = fragment.getSimilarThreshByPos(fracture);
	std::cout << "Segment with fracture:" << fracture << " simThresh: " << simThresh << std::endl;
	std::vector<std::vector<int>> oRegionsList;
	std::vector<std::vector<int>> oRegionOutsideBoundaryVerticesList;
	fragment.segmentByCurvedness(oRegionsList, oRegionOutsideBoundaryVerticesList, simThresh);


	//std::vector<Segment> segments;
	std::vector<Segment> segments; // make it call by reference
	for (int iRegion = 0; iRegion < oRegionsList.size(); iRegion++)
	{
		Segment segment(oRegionsList[iRegion]);
		segment.loadBasicData(fragment.m_VerticesAdjacentFacesList,
			fragment.m_Faces, fragment.m_Vertices,
			fragment.m_Faces2TextureCoordinates, fragment.m_Faces2Normals);
		segments.push_back(segment);
	}


	Visualizer visualizer;
	std::map<int, Eigen::RowVector3d> meshesColors;
	visualizer.generateRandomColors(meshesColors, segments.size());

	for (int i = 0; i < segments.size(); i++)
	{
		//Segment curr = segments[i];
		visualizer.appendMesh(segments[i].m_Vertices, segments[i].m_Faces, meshesColors[i]);
	}

	visualizer.launch();
}