#include "../headers/filters.h"
#include "../headers/lsd.h"
#include "../headers/candidate_discoverer.h"
#include "../headers/cd_params.h"

#include <opencv2/opencv.hpp>

#define AMOUNT_OF_CORNERS 4

double Likelihood(const Candidate& candidate, USE_PARAMS) {
    return PD_K * candidate.perimetr_density
           + AP_K * candidate.area_percentage
           + SP_K * candidate.sides_parallelism;
}

void AddOptionalLines(std::vector<Line>& lines, int32_t width, int32_t height, USE_PARAMS) {
    Line line1{0.0, 0.0,
               static_cast<double>(height) / OPTIONAL_SEGMENT_PERIOD};
    Line line2{0.0, static_cast<double>(width - 1),
               static_cast<double>(height) / OPTIONAL_SEGMENT_PERIOD};
    Line line3{M_PI_2, 0.0,
               static_cast<double>(width) / OPTIONAL_SEGMENT_PERIOD};
    Line line4{M_PI_2, static_cast<double>(height - 1),
               static_cast<double>(width) / OPTIONAL_SEGMENT_PERIOD};
    for (int32_t i = 0; i < height; i += OPTIONAL_SEGMENT_PERIOD * OPTIONAL_SEGMENT_SIZE) {
        line1.segments.emplace_back(Point{i, 0},
                                    Point{i + OPTIONAL_SEGMENT_SIZE, 0});
        line2.segments.emplace_back(Point{i, width - 1},
                                    Point{i + OPTIONAL_SEGMENT_SIZE, width - 1});
    }
    for (int32_t i = 0; i < width; i += OPTIONAL_SEGMENT_PERIOD * OPTIONAL_SEGMENT_SIZE) {
        line3.segments.emplace_back(Point{0, i},
                                    Point{0, i + OPTIONAL_SEGMENT_SIZE});
        line4.segments.emplace_back(Point{height - 1, i},
                                    Point{height - 1, i + OPTIONAL_SEGMENT_SIZE});
    }
    lines.emplace_back(line1);
    lines.emplace_back(line2);
    lines.emplace_back(line3);
    lines.emplace_back(line4);
}

GrayScaledImage GetSingleChanneled(const char* path) {
    cv::Mat img = cv::imread(path, cv::IMREAD_COLOR);
    if (img.empty()) {
        std::cerr << "Error: Could not open or find the image at " << path << std::endl;
        return GrayScaledImage();
    }
    int32_t width = img.cols;
    int32_t height = img.rows;
    Histogram histograms[4];
    for (int y = 0; y < height; y += 2) {
        for (int x = 0; x < width; x += 2) {
            cv::Vec3b color = img.at<cv::Vec3b>(y, x);
            for (int32_t i = 0; i < 3; ++i) {
                histograms[i].Inc(color[i]);
            }
            RGB rgb{color[2], color[1], color[0]};
            histograms[3].Inc(ToGrayScale(rgb));
        }
    }
    int32_t max_var_i = 0;
    double max_var = 0.0;
    for (int32_t i = 0; i < 4; ++i) {
        double var = CalculateHistogramVariance(histograms[i]);
        if (max_var < var) {
            max_var = var;
            max_var_i = i;
        }
    }
    GrayScaledImage single_channeled(width, height);
    for (int y = 0; y < height; ++y) {
        for (int x = 0; x < width; ++x) {
            cv::Vec3b color = img.at<cv::Vec3b>(y, x);
            if (max_var_i < 3) {
                single_channeled.SetPixel(y, x, color[max_var_i]);
            } else {
                RGB rgb{color[2], color[1], color[0]};
                single_channeled.SetPixel(y, x, ToGrayScale(rgb));
            }
        }
    }
    return single_channeled;
}

void ProcessCorners(std::vector<Point>& corners, int32_t reduce_coefficient, int32_t width, int32_t height) {
    for (int32_t i = 0; i < AMOUNT_OF_CORNERS; ++i) {
        corners[i].x *= reduce_coefficient;
        corners[i].y *= reduce_coefficient;
        corners[i].x = std::min(width - 1, std::max(0, corners[i].x));
        corners[i].y = std::min(height - 1, std::max(0, corners[i].y));
    }

    int32_t index_1 = 0;
    double min_dist = ENOUGH_BIG;
    for (int32_t i = 0; i < AMOUNT_OF_CORNERS; ++i) {
        double dist = Distance({0,0}, corners[i]);
        if (dist < min_dist) {
            min_dist = dist;
            index_1 = i;
        }
    }
    std::reverse(corners.begin(), corners.begin() + index_1);
    std::reverse(corners.begin() + index_1, corners.end());
    std::reverse(corners.begin(), corners.end());
    if (CrossProduct(corners[1], corners[0], corners[2]) > 0) {
        std::swap(corners[1], corners[3]);
    }
}

int32_t* DetectCornersImpl(USE_PARAMS) {
    GrayScaledImage single_channeled = GetSingleChanneled(IMG_PATH);
    int32_t width = single_channeled.Width();
    int32_t height = single_channeled.Height();
    int32_t scale_factor = 1;
    GrayScaledImage reduced = Reduce(single_channeled, scale_factor, PARAMS);
    single_channeled.Clear();
    int32_t reduced_width = reduced.Width();
    int32_t reduced_height = reduced.Height();
    GrayScaledImage blured = GaussianBlur(reduced);
    reduced.Clear();
    std::vector<LineSegment> line_segments = DetectLineSegments(blured, PARAMS);
    blured.Clear();
    std::vector<Line> lines = GetLines(line_segments, PARAMS);
    AddOptionalLines(lines, reduced_width, reduced_height, params);
    line_segments.clear();
    auto amount_of_lines = static_cast<int32_t>(lines.size());
    double best_likelihood = 0.0;
    Candidate best_candidate;
    CandidateDiscoverer candidate_discoverer(lines, PARAMS, reduced_width, reduced_height);
    for (int32_t line1 = 0; line1 < amount_of_lines; ++line1) {
        for (int32_t line2 = line1 + 1; line2 < amount_of_lines; ++line2) {
            for (int32_t line3 = line2 + 1; line3 < amount_of_lines; ++line3) {
                for (int32_t line4 = line3 + 1; line4 < amount_of_lines; ++line4) {
                    Candidate candidate = candidate_discoverer.GetCandidate(line1, line2, line3, line4);
                    if (candidate.corners.size() != 4
                        || candidate.area_percentage < 0.1) {
                        continue;
                    }
                    double likelihood = Likelihood(candidate, params);
                    if (likelihood > best_likelihood) {
                        best_likelihood = likelihood;
                        best_candidate = std::move(candidate);
                    }
                }
            }
        }
    }
    ProcessCorners(best_candidate.corners, scale_factor, width, height);
    int32_t* corners = new int32_t[2 * AMOUNT_OF_CORNERS];
    for (int32_t i = 0; i < AMOUNT_OF_CORNERS; ++i) {
        corners[2 * i] = best_candidate.corners[i].x;
        corners[2 * i + 1] = best_candidate.corners[i].y;
    }
    return corners;
}

cv::Point2f Intersection(const cv::Point2f &a1, const cv::Point2f &b1,
                         const cv::Point2f &a2, const cv::Point2f &b2) {
    float x1 = a1.x, y1 = a1.y;
    float x2 = b1.x, y2 = b1.y;
    float x3 = a2.x, y3 = a2.y;
    float x4 = b2.x, y4 = b2.y;
    float denominator = (y4 - y3) * (x2 - x1) - (x4 - x3) * (y2 - y1);
    if (denominator == 0.0f) {
        return {0.0f, 0.0f};
    }
    float ua = ((x4 - x3) * (y1 - y3) - (y4 - y3) * (x1 - x3)) / denominator;
    cv::Point2f pt;
    pt.x = x1 + ua * (x2 - x1);
    pt.y = y1 + ua * (y2 - y1);
    return pt;
}

cv::Point2f ComputeDeviationVector(const std::vector<cv::Point2f>& corners) {
    cv::Point2f center_of_mass(0.0f, 0.0f);
    for (const auto& corner : corners) {
        center_of_mass += corner;
    }
    center_of_mass.x /= static_cast<float>(corners.size());
    center_of_mass.y /= static_cast<float>(corners.size());
    cv::Point2f intersection_point = Intersection(corners[0], corners[2],
                                                  corners[1], corners[3]);
    cv::Point2f deviation_vector = center_of_mass - intersection_point;
    return deviation_vector;
}

cv::Mat TransformDocument(const cv::Mat &input_image, const std::vector<cv::Point2f> &corners, float width, float height) {
    std::vector<cv::Point2f> output_corners = {
            cv::Point2f(0, 0),
            cv::Point2f(width - 1, 0),
            cv::Point2f(width - 1, height - 1),
            cv::Point2f(0, height - 1)
    };
    cv::Mat perspective_matrix = cv::getPerspectiveTransform(corners, output_corners);
    cv::Mat output_image;
    cv::warpPerspective(input_image, output_image, perspective_matrix,
                        cv::Size(static_cast<int32_t>(width), static_cast<int32_t>(height)));
    return output_image;
}

void NormalizeImpl(const char *in_path, const char *out_path, int32_t *corners_ptr) {
    cv::Mat input_image = cv::imread(in_path);
    cv::Point2f c1(static_cast<float>(*(corners_ptr)),
                   static_cast<float>(*(corners_ptr + 1)));
    cv::Point2f c2(static_cast<float>(*(corners_ptr + 2)),
                   static_cast<float>(*(corners_ptr + 3)));
    cv::Point2f c3(static_cast<float>(*(corners_ptr + 4)),
                   static_cast<float>(*(corners_ptr + 5)));
    cv::Point2f c4(static_cast<float>(*(corners_ptr + 6)),
                   static_cast<float>(*(corners_ptr + 7)));

    std::vector <cv::Point2f> corners = {c1, c2, c3, c4};
    cv::Point2f deviation_vector = ComputeDeviationVector(corners);
    float width = static_cast<float>(cv::norm(c1 - c2) + cv::norm(c3 - c4)) / 2
                  + 2 * fabs(deviation_vector.x);
    float height = static_cast<float>(cv::norm(c2 - c3) + cv::norm(c1 - c4)) / 2
                   + 2 * fabs(deviation_vector.y);
    cv::Mat normalized = TransformDocument(input_image, corners, width, height);
    cv::Mat normalized_lab;
    cv::cvtColor(normalized, normalized_lab, cv::COLOR_BGR2Lab);

    std::vector <cv::Mat> lab_channels(3);
    cv::split(normalized_lab, lab_channels);

    cv::Ptr <cv::CLAHE> clahe = cv::createCLAHE();
    clahe->setClipLimit(2.0);
    clahe->setTilesGridSize(cv::Size(8, 8));
    clahe->apply(lab_channels[0], lab_channels[0]);
    cv::Mat normalized_lab_clahe;
    cv::merge(lab_channels, normalized_lab_clahe);
    cv::Mat normalized_bgr_clahe;
    cv::cvtColor(normalized_lab_clahe, normalized_bgr_clahe, cv::COLOR_Lab2BGR);
    cv::imwrite(out_path, normalized_bgr_clahe);
}

extern "C" {
int32_t DetectCorners(const char* path, int32_t* corners_ptr) {
    try {
        INIT_PARAMS
        IMG_PATH = path;
        int32_t* corners_res = DetectCornersImpl(PARAMS);
        for (int32_t i = 0; i < 8; ++i) {
            *(corners_ptr + i) = *(corners_res + i);
        }
        delete[] corners_res;
    } catch(...) {
        return 1;
    }
    return 0;
}

int32_t Normalize(const char* in_path, const char* out_path, int32_t* corners_ptr) {
    try {
        NormalizeImpl(in_path, out_path, corners_ptr);
    } catch(...) {
        return 1;
    }
    return 0;
}
}