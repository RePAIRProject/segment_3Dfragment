#include <string>
#include <vector>
#include "Scripts.h"
#include <Fragment.h>
#include <Segment.h>
#include <string>
#include <iostream>


int main(int argc, char* argv[])
{
	try
	{
		
		std::map<std::string, std::string> params;
		bool isSave = true;
		bool isDebug = false;

		for (int i = 1; i < argc; i++)
		{
			std::string arg = argv[i];

			if (arg == "--disable-save")
			{
				isSave = false;
				continue;
			}

			if (arg == "--enable-debug")
			{
				isDebug = true;
				continue;
			}
			
			if (arg.rfind("--", 0) == 0)
			{
				if (i + 1 < argc)
				{
					std::string value = argv[i + 1];
					params[arg] = value;
					++i;
				}
				else
				{
					std::cerr << "Error: No value provided for parameter " << arg << std::endl;
					return 1;
				}
			}
			else
			{
				std::cerr << "Error: Unknown argument " << arg << std::endl;
				return 1;
			}	
		}
		
		if (params.find("--input-File") == params.end())
		{
			std::cerr << "Error: please provide input file of the fragment to be segmented" << std::endl;
			return 1;
		}

		std::string fragmentPath = params["--input-File"];

		double intactNormalStdThershold = 0.2;
		if (params.find("--intact-Normal-Std-Thershold") != params.end()){ intactNormalStdThershold = std::stod(params["--intact-Normal-Std-Thershold"]); }

		double intactSimilarityFracture = 0.65;
		if (params.find("--intact-Similarity-Fracture") != params.end()) {  intactNormalStdThershold = std::stod(params["--intact-Similarity-Fracture"]); }
		
		std::string script = "segment_intact_surface";
		if (params.find("--script") != params.end()) { script = params["--script"]; }

		
		ObjFragment fragment = ObjFragment(fragmentPath);
		std::cout << "Loading data" << std::endl;
		fragment.load();
		std::cout << "Start script for fragment: " << fragment.m_Name << std::endl;

		if (script == "segment_opposite_surface")
		{
			segment_opposite_surface(fragment, intactNormalStdThershold, intactSimilarityFracture, isSave,isDebug);
		} else
		{
			if (script == "segment_sidewalls_surface")
			{
				segment_sidewalls_surface(fragment, intactNormalStdThershold, intactSimilarityFracture, isSave,isDebug);
			}
			else
			{
				segment_intact_surface(fragment, intactNormalStdThershold, intactSimilarityFracture, isSave, isDebug);
			}
		}
	}
	catch (...)
	{		
		std::cerr << "Caught an unknown exception. A possible workaround: simplify the mesh (using meshlab for example)." << std::endl;
	}
}
