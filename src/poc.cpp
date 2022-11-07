#include <igl/opengl/glfw/Viewer.h>
#include <igl/readOBJ.h>
//#include <igl/png/writePNG.h>
//#include 
//#include "Michal.h"

//#include <igl/barycenter.h>
//#include <igl/eigs.h>

#include <string>
//#include <unistd.h>
#include <Fragment.h>


int old_main(int argc, char* argv[])
{

	//   V  double matrix of vertex positions  #V by 3
 //   TC  double matrix of texture coordinats #TC by 2
 //   N  double matrix of corner normals #N by 3
 //   F  #F list of face indices into vertex positions
 //   FTC  #F list of face indices into vertex texture coordinates
 //   FN  #F list of face indices into vertex normals
	Eigen::MatrixXd V, TC, N;
	Eigen::MatrixXi F, FTC, FN;
	//std::vector<std::tuple<std::string, Index, Index >> FM;

	std::string fragment = "RPf_00154";
	auto bla = igl::readOBJ(".\\fragments\\" + fragment + "\\" + fragment + ".obj", V, TC, N, F, FTC, FN);

	igl::opengl::glfw::Viewer viewer;





	//https://github.com/libigl/libigl/issues/886
	//https://github.com/libigl/libigl/blob/main/include/igl/readOBJ.h
	//https://stackoverflow.com/questions/29867926/why-does-the-number-of-vt-and-v-elements-in-a-blender-obj-file-differ
	//https://stackoverflow.com/questions/27777349/handling-obj-files-why-is-it-possible-to-have-more-vertextextures-vt-than-ve

	auto newVertices = NULL; //std::map<>;

	//std::map<int, double> iglVertices; 
	std::vector<std::vector<double>> iglPositions;
	std::vector<std::vector<double>> iglTextures;
	std::vector<std::vector<double>> iglNormals;
	std::vector<std::vector<double>> colors;
	std::map<std::string, int> uniqueVertices;
	std::vector<std::vector<int>> iglFaces;
	std::vector<int> nextFace;

	for (int faceIndex = 0; faceIndex < F.rows();++faceIndex)
	{

		for (int cornerIndex = 0;cornerIndex < F.cols(); ++cornerIndex)
		{
			int positionIndex = F(faceIndex, cornerIndex);
			int textureIndex = FTC(faceIndex, cornerIndex);
			int normalIndex = FN(faceIndex, cornerIndex);

			std::string nextKey = std::to_string(positionIndex) + "-" + std::to_string(textureIndex) + "-" + std::to_string(normalIndex);

			if (uniqueVertices.count(nextKey) == 0)
			{

				std::vector<double> position = { V(positionIndex,0),V(positionIndex,1) ,V(positionIndex,2) };
				std::vector<double> color = { V(positionIndex,3),V(positionIndex,4) ,V(positionIndex,5) };
				std::vector<double> texture = { TC(textureIndex,0),TC(textureIndex,1) };
				std::vector<double> normal = { N(normalIndex,0),N(normalIndex,1),N(normalIndex,2) };

				iglPositions.push_back(position);
				colors.push_back(color);
				iglTextures.push_back(texture);
				iglNormals.push_back(normal);

				int nextVertIndex = iglPositions.size() - 1;
				uniqueVertices.insert({ nextKey,nextVertIndex });

				nextFace.push_back(nextVertIndex);
			}
			else {
				nextFace.push_back(uniqueVertices[nextKey]);
			}
		}
		iglFaces.push_back(nextFace);
		nextFace.clear();
	}

	int numUniqueComb = iglPositions.size();
	Eigen::MatrixXd iglPositionsMatrix = Eigen::MatrixXd::Zero(numUniqueComb, 3);
	Eigen::MatrixXd iglColorsMatrix = Eigen::MatrixXd::Zero(numUniqueComb, 3);
	Eigen::MatrixXd iglTextureMatrix = Eigen::MatrixXd::Zero(numUniqueComb, 2);
	Eigen::MatrixXd iglNormalsMatrix = Eigen::MatrixXd::Zero(numUniqueComb, 3);

	for (int combIndex = 0; combIndex < numUniqueComb;++combIndex)
	{

		auto pos = iglPositions.begin() + combIndex;
		auto color = colors.begin() + combIndex;
		auto tex = iglTextures.begin() + combIndex;
		auto norm = iglNormals.begin() + combIndex;
		int i = 0;

		for (auto coord = pos->begin(); coord != pos->end(); ++coord)
		{
			iglPositionsMatrix(combIndex, i) = *coord;
			++i;
		}

		i = 0;
		for (auto coord = color->begin(); coord != color->end(); ++coord)
		{
			iglColorsMatrix(combIndex, i) = *coord;
			++i;
		}

		i = 0;
		for (auto coord = tex->begin(); coord != tex->end(); ++coord)
		{
			iglTextureMatrix(combIndex, i) = *coord;
			++i;
		}

		i = 0;
		for (auto coord = norm->begin(); coord != norm->end(); ++coord)
		{
			iglNormalsMatrix(combIndex, i) = *coord;
			++i;
		}
	}

	int numFaces = iglFaces.size();
	Eigen::MatrixXi iglFacesMatrix = Eigen::MatrixXi::Zero(numFaces, 3);

	int i = 0;
	for (int faceIndex = 0; faceIndex < numFaces; ++faceIndex)
	{
		auto face = iglFaces.begin() + faceIndex;
		i = 0;
		for (auto vertIndex = face->begin(); vertIndex != face->end(); ++vertIndex)
		{
			iglFacesMatrix(faceIndex, i) = *vertIndex;
			++i;
		}
	}




	/// try to read ply object
	Eigen::MatrixXd ply_v, ply_colors, ply_N, ply_UV, ply_vd, ply_FD, ply_ED;
	Eigen::MatrixXi ply_f, ply_E;
	std::vector<std::string> ply_Vheader, ply_Eheader, ply_Fheader, ply_comments;
	igl::readPLY(".\\fragments\\" + fragment + "\\" + fragment + ".ply", ply_v, ply_f, ply_E, ply_N, ply_UV,
		ply_vd, ply_Vheader, ply_FD, ply_Fheader, ply_ED, ply_Eheader, ply_comments);

	if (ply_vd.cols() > 1)
	{
		Eigen::MatrixXd sub_color_matrix = ply_vd.block(0, 0, ply_vd.rows(), 3);
		if (sub_color_matrix.size() != 0)
		{
			ply_colors =
				(sub_color_matrix.rowwise() - sub_color_matrix.colwise().minCoeff()).array().rowwise() /
				(sub_color_matrix.colwise().maxCoeff() - sub_color_matrix.colwise().minCoeff()).array();
		}
	}
	else //no texture on object
	{
		ply_colors = Eigen::MatrixXd::Zero(ply_v.rows(), ply_v.cols());
	}

	//viewer.append_mesh(false);
	////viewer.data().set_vertices(v);
	//viewer.data_list[0].set_mesh(V.block(0, 0, V.rows(), 3), F);
	//viewer.data().set_face_based(true);
	//viewer.data_list[0].set_colors(V.block(0, 3, V.rows(), 3));
	//viewer.data_list[0].set_normals(N);
	////viewer.data().set_uv(TC, FTC);
	//viewer.data_list[0].show_lines = false;
	////viewer.data_list[0].show_texture = true;

	////viewer.append_mesh();//viewer.append_mesh(false);
	viewer.data_list[0].set_mesh(iglPositionsMatrix, iglFacesMatrix);
	viewer.data_list[0].set_normals(iglNormalsMatrix);
	viewer.data_list[0].set_uv(iglTextureMatrix);
	viewer.data_list[0].set_colors(iglColorsMatrix);
	viewer.data_list[0].show_lines = false;
	viewer.data_list[0].show_texture = true;
	//viewer.data_list[0].use_matcap = true;
	/*viewer.data_list[0].uniform_colors(Eigen::Vector3d(0.200000,0.200000,0.200000),
										Eigen::Vector3d(1.000000,1.000000,1.000000),
										Eigen::Vector3d(1.000000, 1.000000, 1.000000));*/


										/*//viewer.append_mesh(false);
										viewer.data_list[0].set_mesh(ply_v, ply_f);
										viewer.data_list[0].set_colors(ply_colors);
										viewer.data_list[0].show_lines = false;*/
										////viewer.data_list[2].show_texture = true;




	viewer.callback_key_down = [&](igl::opengl::glfw::Viewer&, unsigned int key, int mod) -> bool {
		if (key == GLFW_KEY_L)
		{
			viewer.data_list[viewer.selected_data_index].set_visible(false);
			viewer.selected_data_index = (viewer.selected_data_index + 1) % viewer.data_list.size();
			viewer.data_list[viewer.selected_data_index].set_visible(true);
			viewer.data_list[viewer.selected_data_index].show_lines = false;
			std::cout << "The selected is " << viewer.selected_data_index << std::endl;
		}

		if (key == GLFW_KEY_J)
		{
			viewer.core().lighting_factor = 2 * viewer.core().lighting_factor;
		}

		if (key == GLFW_KEY_K)
		{
			viewer.core().lighting_factor = 0.5 * viewer.core().lighting_factor;
		}

		if (key == GLFW_KEY_Z)
		{
			if (viewer.core().lighting_factor == 0)
			{
				viewer.core().lighting_factor = 1;
			}
			else
			{
				viewer.core().lighting_factor = 0;

			}

		}

		if (key == GLFW_KEY_S)
		{
			igl::writeOBJ(".\\fragments\\" + fragment + "\\" + fragment + "_igl.obj", V, F, N, FN, TC, FTC);
		}

		if (key == GLFW_KEY_C)
		{
			Eigen::MatrixXd pd1, pd2;
			Eigen::MatrixXd pv1, pv2;

			//igl::principal_curvature(V, F, pd1, pd2, pv1, pv2);
		}

		return true;
	};




	viewer.launch();
	return 0;
}
