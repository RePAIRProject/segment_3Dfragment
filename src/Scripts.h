#pragma once
#include <string>
#include <vector>
#include <Fragment.h>
#include <Segment.h>


Segment segment_intact_surface(ObjFragment &fragment, std::string outFileName);
void segment_opposite_surface(ObjFragment& fragment);