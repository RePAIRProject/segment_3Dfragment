#pragma once
#include <string>
#include <vector>
#include <Fragment.h>
#include <Segment.h>
#include <Eigen/Eigenvalues>

Segment segmentIntactSurface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave, bool isVisualizer);
void segmentOppositeSurface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave,bool isVisualizer);
void segmentSidewallsSurface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave, bool isVisualizer);