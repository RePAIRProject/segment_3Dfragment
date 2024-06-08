#pragma once
#include <Segment.h>
#include <Fragment.h>
#include <igl/writeOBJ.h>
#include <iostream>
#include <cstdio>
#include <filesystem>



void initSegments(std::vector<Segment> &oSegments, std::map<int, int> &oVertIndex2SegIndex,
					const std::vector<std::vector<int>> &oRegionsList,
					const std::vector<std::vector<int>> &oRegionOutsideBoundaryVerticesList,
					 double fragmentSize, ObjFragment& parentFragment);

void sortToSmallAndBigSegments(std::map<int, Segment*>& oSmallSegments, std::map<int, Segment*>& oBigSegments, 
							std::vector<Segment> &segments, double minBigSegPercSize);

Eigen::Vector3d calcAvg(const std::vector<Eigen::Vector3d>& vectors);
Eigen::Vector3d calcVariance(const std::vector<Eigen::Vector3d>& vectors, Eigen::Vector3d mean);
void merge(int iSrcSeg, int iDstSeg, std::map<int, int>& vertIndex2SegIndex, std::map<int, Segment*>& srcSegPool, std::map<int, Segment*>& dstSegPool);
void merge(Segment* srcSeg, Segment* dstSeg, int iDstSeg, std::map<int, int>& vertIndex2SegIndex);
void mergeSmall2BigSegments(std::map<int, Segment*>& smallSegments, std::map<int, Segment*>& bigSegments, std::map<int, int>& vertIndex2SegIndex);

double sigmiod(double x);


void saveMtlFile(std::string outPath, std::string img, std::string materialName="material_0");

void saveObjFile(std::string outObjPath, std::string outMtlPath,
	const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, const Eigen::MatrixXd& normals,
	const Eigen::MatrixXi& faces2Normals, const Eigen::MatrixXd& textureCoordinates,
	const Eigen::MatrixXi& faces2TextureCoordinates, std::string materialName = "material_0");