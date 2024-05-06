#include "../headers/candidate_discoverer.h"

#include <algorithm>

Candidate::Candidate() : corners(4) {
    for (int32_t i = 0; i < 4; ++i) {
        corners[i] = Point{ENOUGH_BIG, ENOUGH_BIG};
    }
}

CandidateDiscoverer::CandidateDiscoverer(const std::vector<Line> &lines, USE_PARAMS,
                                         int32_t img_width, int32_t img_height)
        : line_intersections_(static_cast<int32_t>(lines.size()), -1),
          line_angles_(lines.size()), img_width_(img_width), img_height_(img_height),
          prefix_weights_(lines.size()) {
    for (int32_t i = 0; i < lines.size(); ++i) {
        line_angles_[i] = lines[i].angle;
    }
    std::vector<std::vector<int32_t>> points_on_lines(lines.size());
    for (int32_t i = 0; i < lines.size(); ++i) {
        for (int32_t j = i + 1; j < lines.size(); ++j) {
            Point intersection = Intersect(lines[i], lines[j]);
            if (intersection.x < -AIRBAG || intersection.x > img_width + AIRBAG
                || intersection.y < -AIRBAG || intersection.y > img_height + AIRBAG) {
                continue;
            }
            points_on_lines[i].emplace_back(points_.size());
            points_on_lines[j].emplace_back(points_.size());
            line_intersections_[i][j] = static_cast<int32_t>(points_.size());
            line_intersections_[j][i] = static_cast<int32_t>(points_.size());
            points_.emplace_back(intersection);
        }
    }
    ComparatorParam param;
    auto cmp = [this, &param] (int32_t p1_index, int32_t p2_index) {
        if (param == ComparatorParam::ByX) {
            return points_[p1_index].x < points_[p2_index].x;
        }
        return points_[p1_index].y < points_[p2_index].y;
    };
    for (int32_t i = 0; i < lines.size(); ++i) {
        if (i == 3) {
            int z = 1 + 1;
        }
        double s = sin(lines[i].angle);
        double c = cos(lines[i].angle);
        if (fabs(c) < fabs(s)) {
            param = ComparatorParam::ByX;
        } else {
            param = ComparatorParam::ByY;
        }
        std::sort(points_on_lines[i].begin(), points_on_lines[i].end(), cmp);

        int32_t seg_index = 0;
        double total_length = 0.0;
        for (auto point_index : points_on_lines[i]) {
            Point point = points_[point_index];
            while (seg_index < lines[i].segments.size() &&
                   ((param == ByX) ? (lines[i].segments[seg_index].second.x < point.x)
                                   : (lines[i].segments[seg_index].second.y < point.y))) {
                total_length += Distance(lines[i].segments[seg_index].first,
                                         lines[i].segments[seg_index].second);
                ++seg_index;
            }
            if (seg_index < lines[i].segments.size() &&
                ((param ==ByX) ? (lines[i].segments[seg_index].first.x < point.x)
                               : (lines[i].segments[seg_index].first.y < point.y) )) {
                prefix_weights_[i][point_index] = total_length
                                                  + Distance(point, lines[i].segments[seg_index].first);
            } else {
                prefix_weights_[i][point_index] = total_length;
            }
        }
    }
}

Candidate CandidateDiscoverer::GetCandidate(int32_t line1, int32_t line2, int32_t line3, int32_t line4) {
    Candidate candidate;
    int32_t i1, i2, i3, i4;
    for (int32_t i = 0; i < 3; ++i) {
        std::swap(line1, line2);
        int32_t l = line1;
        line1 = line2;
        line2 = line3;
        line3 = line4;
        line4 = l;
        i1 = line_intersections_[line1][line2];
        i2 = line_intersections_[line2][line3];
        i3 = line_intersections_[line3][line4];
        i4 = line_intersections_[line4][line1];
        if (i1 == -1 || i2 == -1 || i3 == -1 || i4 == -1) {
            continue;
        }
        candidate.corners[0] = points_[i1];
        candidate.corners[1] = points_[i2];
        candidate.corners[2] = points_[i3];
        candidate.corners[3] = points_[i4];
        if (IsConvex(candidate.corners[0], candidate.corners[1],
                     candidate.corners[2], candidate.corners[3])
            && IsNotSelfIntersected(candidate.corners[0], candidate.corners[1],
                                    candidate.corners[2], candidate.corners[3])) {
            break;
        }
    }
    if (candidate.corners[0].x == ENOUGH_BIG || !IsConvex(candidate.corners[0], candidate.corners[1],
                  candidate.corners[2], candidate.corners[3])
                  || !IsNotSelfIntersected(candidate.corners[0], candidate.corners[1],
                                           candidate.corners[2], candidate.corners[3])) {
        candidate.corners.clear();
        return candidate;
    }
    candidate.sides_parallelism = (fabs(cos(line_angles_[line1] - line_angles_[line3]))
                                    + fabs(cos(line_angles_[line2] - line_angles_[line4]))) / 2;
    double area = GetArea(candidate.corners);
    candidate.area_percentage = area / (img_height_ * img_width_);
    candidate.perimetr_density += fabs(prefix_weights_[line2][i2] - prefix_weights_[line2][i1]);
    candidate.perimetr_density += fabs(prefix_weights_[line3][i3] - prefix_weights_[line3][i2]);
    candidate.perimetr_density += fabs(prefix_weights_[line4][i4] - prefix_weights_[line4][i3]);
    candidate.perimetr_density += fabs(prefix_weights_[line1][i1] - prefix_weights_[line1][i4]);
    double perimetr = Distance(candidate.corners[0], candidate.corners[1])
                    + Distance(candidate.corners[1], candidate.corners[2])
                    + Distance(candidate.corners[2], candidate.corners[3])
                    + Distance(candidate.corners[3], candidate.corners[0]);
    candidate.perimetr_density /= perimetr;
    return candidate;
}
