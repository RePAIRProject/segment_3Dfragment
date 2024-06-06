#pragma once
#include <string>
#include <vector>
#include <Fragment.h>
#include <Segment.h>
#include <Eigen/Eigenvalues>

Segment segment_intact_surface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave, bool isVisualizer);
void segment_opposite_surface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave,bool isVisualizer);
void segment_sidewalls_surface(ObjFragment& fragment, double intactNormalStdThershold, double intactSimilarityFracture, bool isSave, bool isVisualizer);