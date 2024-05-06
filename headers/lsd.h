#pragma once

#include "image.h"
#include "helpers.h"
#include "cd_params.h"

std::vector<LineSegment> DetectLineSegments(const GrayScaledImage& img, USE_PARAMS);
std::vector<Line> GetLines(std::vector<LineSegment>& line_segments, USE_PARAMS);
