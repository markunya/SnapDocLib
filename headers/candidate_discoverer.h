#pragma once

#include "helpers.h"
#include "cd_params.h"

#include <unordered_map>

struct Candidate {
    std::vector<Point> corners;
    double perimetr_density = 0.0;
    double area_percentage = 0.0;
    double sides_parallelism = 0.0;
    Candidate();
};

class CandidateDiscoverer {
public:
    CandidateDiscoverer(const std::vector<Line>& lines, USE_PARAMS,
                        int32_t img_width, int32_t img_height);

    Candidate GetCandidate(int32_t line1, int32_t line2, int32_t line3, int32_t line4);

private:
    Matrix line_intersections_;
    std::vector<double> line_angles_;
    int32_t img_width_;
    int32_t img_height_;
    std::vector<Point> points_;
    std::vector<std::unordered_map<int32_t, double>> prefix_weights_;
};
