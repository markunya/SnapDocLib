#include "../headers/lsd.h"

#include "../headers/helpers.h"

#include <cmath>
#include <queue>
#include <algorithm>

#define BIN_SIZE 1024
#define LogGamma(x) ((x) > 15.0 ? LogGammaWindshcitl(x) : LogGammaLanczos(x))
#define TABSIZE 100000

static double LogGammaLanczos(double x) {
  static double q[7] = { 75122.6331530, 80916.6278952, 36308.2951477,
                         8687.24529705, 1168.92649479, 83.8676043424,
                         2.50662827511 };
  double a = (x + 0.5) * log(x + 5.5) - (x + 5.5);
  double b = 0.0;
  for(int32_t i = 0; i < 7; ++i) {
      a -= log(x + static_cast<double>(i));
      b += q[i] * pow(x, static_cast<double>(i));
    }
  return a + log(b);
}

static double LogGammaWindshcitl(double x) {
  return 0.918938533204673 + (x - 0.5) * log(x) - x
         + 0.5 * x * log(x * sinh(1 / x) + 1 / (810.0 * pow(x,6.0)));
}

static bool IsAligned(double a, double b, double tolerance) {
    while (a < -M_PI) {
        a += M2_PI;
    }
    while (a > M2_PI) {
        a -= M2_PI;
    }
    while (b < -M_PI) {
        b += M2_PI;
    }
    while (b > M_PI) {
        b -= M2_PI;
    }
    if (fabs(a - b) < tolerance) {
        return true;
    }
    if (fabs(fabs(a - b) - M2_PI) < tolerance) {
        return true;
    }
    return false;
}

static double SignedAngleDifference(double a, double b) {
    a -= b;
    while( a <= -M_PI ) {
        a += M2_PI;
    }
    while( a > M_PI ) {
        a -= M2_PI;
    }
    return a;
}

static double AngleDifference(double a, double b)
{
    a -= b;
    while( a <= -M_PI ) {
        a += M2_PI;
    }
    while( a > M_PI ) {
        a -= M2_PI;
    }
    if( a < 0.0 ) a = -a;
    return a;
}

std::vector<std::vector<Point>> GetLineSegmentAngles(const GrayScaledImage& img, BoolImage& used,
                          DoubleImage& gradient_angles, DoubleImage& gradient_magnitudes,
                          double gradient_magnitude_threshold) {
    std::vector<std::vector<Point>> ordered_indexes(BIN_SIZE);
    double max_magnitude = 0.0;
    for (int32_t i = 0; i < img.Height(); ++i) {
        for (int32_t j = 0; j < img.Width(); ++j) {
            int32_t dx = (img.GetPixel(i, j + 1) + img.GetPixel(i + 1, j + 1)
                          - img.GetPixel(i, j) - img.GetPixel(i + 1, j));
            int32_t dy = (img.GetPixel(i + 1, j) + img.GetPixel(i + 1, j + 1)
                          - img.GetPixel(i, j) - img.GetPixel(i, j + 1));
            double magnitude = sqrt((dx * dx + dy * dy) / 4.0);
            gradient_magnitudes.SetPixel(i, j, magnitude);
            if (magnitude > max_magnitude) {
                max_magnitude = magnitude;
            }
            if (magnitude <= gradient_magnitude_threshold) {
                used.SetPixel(i, j, true);
                continue;
            }
            gradient_angles.SetPixel(i, j, atan2(dx, -dy));
        }
    }
    for (int32_t i = 0; i < gradient_magnitudes.Height(); ++i) {
        for (int32_t j = 0; j < gradient_magnitudes.Width(); ++j) {
            auto index = static_cast<int32_t>(gradient_magnitudes.GetPixel(i, j) * BIN_SIZE / max_magnitude);
            if (index >= BIN_SIZE) {
                index = BIN_SIZE - 1;
            }
            Point p = {i, j};
            ordered_indexes[BIN_SIZE - 1 - index].emplace_back(p);
        }
    }
    return ordered_indexes;
}

struct Region {
    double angle = 0.0;
    std::vector<Point> points;
};

static Region GrowRegion(const DoubleImage& angles,
                         BoolImage& used,
                         double angle_tolerance,
                         Point start_point) {
    Region region;
    int32_t width = angles.Width();
    int32_t height = angles.Height();
    int32_t dx[8] = {-1, 0 , 1, 1, 1, 0, -1, -1};
    int32_t dy[8] = {1, 1, 1, 0, -1, -1, -1, 0};
    auto [start_y, start_x] = start_point;
    region.angle = angles.GetPixel(start_y, start_x);
    double sum_sin = 0.0;
    double sum_cos = 0.0;
    std::queue<Point> bfs_queue;
    bfs_queue.push(Point{start_y, start_x});
    while (!bfs_queue.empty()) {
        auto [y, x] = bfs_queue.front();
        bfs_queue.pop();
        if (y < 0 || x < 0 || y >= height || x >= width) {
            continue;
        }
        if (used.GetPixel(y, x)) {
            continue;
        }
        used.SetPixel(y, x, true);
        if (!IsAligned(region.angle, angles.GetPixel(y, x), angle_tolerance)) {
            continue;
        }
        sum_sin += sin(angles.GetPixel(y, x));
        sum_cos += cos(angles.GetPixel(y, x));
        region.angle = atan2(sum_sin, sum_cos);
        region.points.emplace_back(Point{y, x});
        for (int32_t k = 0; k < 8; ++k) {
            bfs_queue.emplace(Point{y + dy[k], x + dx[k]});
        }
    }
    return region;
}

static Rectangle RegionToRectangle(const Region& region, const DoubleImage& magnitudes, double angle_tolerance, double probability) {
    Rectangle rectangle;
    double sum_of_magnitudes = 0.0;
    for (auto [y, x] : region.points) {
        sum_of_magnitudes += magnitudes.GetPixel(y, x);
        rectangle.center_x += magnitudes.GetPixel(y, x) * static_cast<double>(x);
        rectangle.center_y += magnitudes.GetPixel(y, x) * static_cast<double>(y);
    }
    rectangle.center_x /= sum_of_magnitudes;
    rectangle.center_y /= sum_of_magnitudes;
    double Mxx = 0.0;
    double Myy = 0.0;
    double Mxy = 0.0;
    for (auto [y, x] : region.points) {
        Mxx += (rectangle.center_x - static_cast<double>(x)) * (rectangle.center_x - static_cast<double>(x)) * magnitudes.GetPixel(y, x);
        Myy += (rectangle.center_y - static_cast<double>(y)) * (rectangle.center_x - static_cast<double>(y)) * magnitudes.GetPixel(y, x);
        Mxy += (rectangle.center_y - static_cast<double>(y)) * (rectangle.center_x - static_cast<double>(x)) * magnitudes.GetPixel(y, x);
    }
    Mxx /= sum_of_magnitudes;
    Myy /= sum_of_magnitudes;
    Mxy /= sum_of_magnitudes;
    double lambda = (Mxx + Myy + sqrt((Mxx - Myy) * (Mxx - Myy) + 4.0 * Mxy * Mxy)) / 2.0;
    rectangle.angle = (fabs(Mxx) > fabs(Myy) ? atan2(lambda - Mxx,Mxy) : atan2(Mxy,lambda - Myy));
    if(AngleDifference(rectangle.angle, region.angle) > angle_tolerance ) {
        if (rectangle.angle > 0) {
            rectangle.angle -= M_PI;
        } else {
            rectangle.angle += M_PI;
        }
    }
    rectangle.vx = cos(rectangle.angle);
    rectangle.vy = sin(rectangle.angle);
    double l_min = 0.0;
    double l_max = 0.0;
    double w_min = 0.0;
    double w_max = 0.0;
    for(int32_t  i =0; i < region.points.size(); ++i) {
        double l = (static_cast<double>(region.points[i].x) - rectangle.center_x) * rectangle.vx
                + (static_cast<double>(region.points[i].y) - rectangle.center_y) * rectangle.vy;
        double w = -(static_cast<double>(region.points[i].x) - rectangle.center_x) * rectangle.vy
                + (static_cast<double>(region.points[i].y) - rectangle.center_y) * rectangle.vx;
        if( l > l_max ) {
            l_max = l;
        }
        if( l < l_min ) {
            l_min = l;
        }
        if( w > w_max ) {
            w_max = w;
        }
        if( w < w_min ) {
            w_min = w;
        }
    }
    rectangle.x1 = rectangle.center_x + l_min * rectangle.vx;
    rectangle.y1 = rectangle.center_y + l_min * rectangle.vy;
    rectangle.x2 = rectangle.center_x + l_max * rectangle.vx;
    rectangle.y2 = rectangle.center_y + l_max * rectangle.vy;
    rectangle.width = (w_max - w_min);
    rectangle.angle_tolerance = angle_tolerance;
    rectangle.probability = probability;
    return rectangle;
}

static bool ReduceRegionRadius(Region& region, const DoubleImage& magnitudes,
                                 double angle_tolerance, Rectangle& rectangle,
                                 BoolImage& used, double req_density, double probability,
                                 int32_t min_reg_size) {
    double density = static_cast<double>(region.points.size()) /
              (Distance(rectangle.x1,rectangle.y1,rectangle.x2,rectangle.y2) * rectangle.width);
    if(density >= req_density){
        return true;
    }

    auto x = static_cast<double>(region.points[0].x);
    auto y = static_cast<double>(region.points[0].y);
    for (auto& [ya, xa] : region.points) {
        used.SetPixel(ya, xa, false);
    }
    std::sort(region.points.begin(), region.points.end(), [x, y](const auto& p1, const auto& p2) {
        return Distance(x, y, p1.x, p1.y) < Distance(x, y, p2.x, p2.y) ;
    });
    while(density < req_density) {
        region.points.resize(static_cast<size_t>(static_cast<double>(region.points.size()) * 0.75));

        if(region.points.size() < min_reg_size){
            for (auto& [ya, xa] : region.points) {
                used.SetPixel(ya, xa, true);
            }
            return false;
        }

        rectangle = RegionToRectangle(region, magnitudes, angle_tolerance, probability);

        density = static_cast<double>(region.points.size()) /
                  (Distance(rectangle.x1, rectangle.y1, rectangle.x2, rectangle.y2) * rectangle.width);
    }
    for (auto& [ya, xa] : region.points) {
        used.SetPixel(ya, xa, true);
    }
    return true;
}

static bool RefineRectangle(const Rectangle& rectangle, const Region& region,
                     const DoubleImage& magnitudes, const DoubleImage& angles, BoolImage& used,
                     double angle_tolerance, double req_density, double probability, int32_t min_reg_size) {
    double density = static_cast<double>(region.points.size()) /
                     (rectangle.width * Distance(rectangle.x1, rectangle.y1, rectangle.x2, rectangle.y2));
    if (density >= req_density) {
        return true;
    }
    auto [y, x] = region.points[0];
    double angle = angles.GetPixel(y, x);
    double sum_of_diffs = 0.0;
    double sum_of_squared_diffs = 0.0;
    int32_t amount = 0;
    for (auto [ya, xa]: region.points) {
        used.SetPixel(ya, xa, false);
        if (Distance(x, y, xa, ya) < rectangle.width) {
            double angle_a = angles.GetPixel(ya, xa);
            double angle_diff = SignedAngleDifference(angle_a, angle);
            sum_of_diffs += angle_diff;
            sum_of_squared_diffs += angle_diff * angle_diff;
            ++amount;
        }
    }
    double mean_angle = sum_of_diffs / static_cast<double>(amount);
    double tau = 2.0 * sqrt((sum_of_squared_diffs - 2.0 * mean_angle * sum_of_diffs)
                            / static_cast<double>(amount) + mean_angle * mean_angle);

    Region region2 = GrowRegion(angles, used, tau, region.points[0]);

    if (region2.points.size() < 2) {
        return false;
    }

    Rectangle rectangle2 = RegionToRectangle(region2, magnitudes, angle_tolerance, probability);

    density = static_cast<double>(region2.points.size()) /
              (Distance(rectangle2.x1, rectangle.y1, rectangle.x2, rectangle.y2) * rectangle.width);

    if (density < req_density) {
        return ReduceRegionRadius(region2, magnitudes, angle_tolerance,rectangle2,
                                  used, req_density, probability, min_reg_size);
    }
    return true;
}

static double NFA(int32_t amount, int32_t amount_of_aligned, double probability, double log_number_of_tests) {
  static double inv[TABSIZE];
  double tolerance = 0.1;
  double log1term, term, bin_term, mult_term, bin_tail, err, p_term;

  if(amount == 0 || amount_of_aligned == 0){
      return -log_number_of_tests;
  }
  if(amount_of_aligned == amount){
      return -log_number_of_tests - static_cast<double>(amount) * log10(probability);
  }

  p_term = probability / (1.0 - probability);

  log1term = LogGamma(static_cast<double>(amount) + 1.0) - LogGamma(static_cast<double>(amount_of_aligned) + 1.0)
           - LogGamma(static_cast<double>(amount - amount_of_aligned) + 1.0 )
           + static_cast<double>(amount_of_aligned) * log(probability)
           + static_cast<double>(amount - amount_of_aligned) * log(1.0 - probability);
  term = exp(log1term);

  if(DoubleEqual(term,0.0)) {
      if(static_cast<double>(amount_of_aligned) > static_cast<double>(amount) * probability) {
          return -log1term / M_LN10 - log_number_of_tests;
      }
      else
        return -log_number_of_tests;
    }

  bin_tail = term;
  for(int32_t i = amount_of_aligned + 1; i <= amount; ++i) {
      bin_term = static_cast<double>(amount - i + 1) * (i < TABSIZE ?
                   (inv[i] != 0.0 ? inv[i] : (inv[i] = 1.0 / static_cast<double>(i))) :
                   1.0 / static_cast<double>(i));

      mult_term = bin_term * p_term;
      term *= mult_term;
      bin_tail += term;
      if(bin_term < 1.0) {
          err = term * (( 1.0 - pow( mult_term, static_cast<double>(amount - i + 1))) /
                         (1.0-mult_term) - 1.0);
          if(err < tolerance * fabs(-log10(bin_tail) - log_number_of_tests) * bin_tail) {
              break;
          }
        }
    }
  return -log10(bin_tail) - log_number_of_tests;
}

static double RectangleNFA(const Rectangle& rectangle, const DoubleImage& angles, double log_number_of_tests) {
    RectangleIterator it(rectangle);
    int32_t amount = 0;
    int32_t amount_of_aligned = 0;
    while (!it.IsEnd()) {
        auto [y, x] = it.GetPoint();
        ++amount;
        if (IsAligned(rectangle.angle, angles.GetPixel(y, x), rectangle.angle_tolerance)) {
            ++amount_of_aligned;
        }
        it.Next();
    }
    return NFA(amount, amount_of_aligned, rectangle.probability, log_number_of_tests);
}

static double ImproveRectangle(Rectangle& rectangle, const DoubleImage& angles,
                        double log_number_of_tests, double log_eps) {
   Rectangle r;
   double log_nfa,log_nfa_new;
   double delta = 0.5;
   double delta_2 = delta / 2.0;

   log_nfa = RectangleNFA(rectangle, angles, log_number_of_tests);
   if(log_nfa > log_eps){
       return log_nfa;
   }

  r = rectangle;
  for(int32_t i = 0; i < 5; ++i) {
      r.probability /= 2.0;
      r.angle_tolerance = r.probability * M_PI;
      log_nfa_new = RectangleNFA(r, angles, log_number_of_tests);
      if(log_nfa_new > log_nfa) {
          log_nfa = log_nfa_new;
          rectangle = r;
        }
    }

  if(log_nfa > log_eps) {
      return log_nfa;
  }

  r = rectangle;
  for(int32_t i = 0; i < 5; ++i) {
      if((r.width - delta) >= 0.5) {
          r.width -= delta;
          log_nfa_new = RectangleNFA(r, angles, log_number_of_tests);
          if(log_nfa_new > log_nfa) {
              rectangle = r;
              log_nfa = log_nfa_new;
          }
      }
  }

  if( log_nfa > log_eps ) {
      return log_nfa;
  }

  r = rectangle;
  for(int32_t i = 0; i < 5; ++i) {
      if((r.width - delta) >= 0.5) {
          r.x1 += -r.vy * delta_2;
          r.y1 +=  r.vx * delta_2;
          r.x2 += -r.vy * delta_2;
          r.y2 +=  r.vx * delta_2;
          r.width -= delta;
          log_nfa_new = RectangleNFA(r,angles, log_number_of_tests);
          if(log_nfa_new > log_nfa) {
              rectangle = r;
              log_nfa = log_nfa_new;
          }
      }
  }

  if(log_nfa > log_eps){
      return log_nfa;
  }

  r = rectangle;
  for(int32_t i = 0; i < 5; ++i) {
      if((r.width - delta) >= 0.5) {
          r.x1 -= -r.vy * delta_2;
          r.y1 -=  r.vx * delta_2;
          r.x2 -= -r.vy * delta_2;
          r.y2 -=  r.vx * delta_2;
          r.width -= delta;
          log_nfa_new = RectangleNFA(r, angles,log_number_of_tests);
          if(log_nfa_new > log_nfa) {
              rectangle = r;
              log_nfa = log_nfa_new;
          }
      }
  }

  if( log_nfa > log_eps ){
      return log_nfa;
  }

  r = rectangle;
  for(int32_t i = 0; i < 5; ++i) {
      r.probability /= 2.0;
      r.angle_tolerance = r.probability * M_PI;
      log_nfa_new = RectangleNFA(r, angles, log_number_of_tests);
      if(log_nfa_new > log_nfa) {
          log_nfa = log_nfa_new;
          rectangle = r;
      }
  }

  return log_nfa;
}

std::vector<LineSegment> DetectLineSegments(const GrayScaledImage& img, USE_PARAMS) {
    std::vector<LineSegment> result;
    double log_number_of_tests = 5.0 * ( log10( img.Width()) + log10( img.Height()) ) / 2.0
                                 + log10(11.0);
    auto min_reg_size = MIN_REG_SIZE_K * static_cast<int32_t>(-log_number_of_tests
                                                              / log10(ANGLE_TOLERANCE / M_PI));
    BoolImage used(img.Width(), img.Height());
    DoubleImage gradient_angles(img.Width(), img.Height());
    DoubleImage gradient_magnitudes(img.Width(), img.Height());
    auto ordered_indexes =
            GetLineSegmentAngles(img, used, gradient_angles,
                                 gradient_magnitudes, GRAD_MAGNITUDE_THRESHOLD);
    for (int32_t i = 0; i < BIN_SIZE; ++i) {
        for (auto [y, x] : ordered_indexes[i]) {
            if (used.GetPixel(y, x)) {
                continue;
            }
            Region region = GrowRegion(gradient_angles, used,
                                       ANGLE_TOLERANCE,
                                       Point{y, x});
            if (region.points.size() < min_reg_size) {
                continue;
            }
            Rectangle rectangle = RegionToRectangle(region, gradient_magnitudes,
                                                    ANGLE_TOLERANCE, PROBABILITY);
            if (!RefineRectangle(rectangle, region, gradient_magnitudes,gradient_angles,
                                 used, ANGLE_TOLERANCE, REQ_DENSITY,
                                 PROBABILITY, min_reg_size)) {
                continue;
            }
            double log_nfa = ImproveRectangle(rectangle, gradient_angles,
                                              log_number_of_tests, LOG_EPS);
            if(log_nfa <= LOG_EPS) {
                continue;
            }
            rectangle.x1 += 0.5;
            rectangle.y1 += 0.5;
            rectangle.x2 += 0.5;
            rectangle.y2 += 0.5;
            LineSegment line_segment;
            line_segment.x1 = std::min(img.Width() - 1, std::max(0, static_cast<int32_t>(rectangle.x1)));
            line_segment.y1 = std::min(img.Height() - 1, std::max(0, static_cast<int32_t>(rectangle.y1)));
            line_segment.x2 = std::min(img.Width() - 1, std::max(0, static_cast<int32_t>(rectangle.x2)));
            line_segment.y2 = std::min(img.Height() - 1, std::max(0, static_cast<int32_t>(rectangle.y2)));
            line_segment.angle = rectangle.angle;
            line_segment.length = Distance(line_segment.x1, line_segment.y1, line_segment.x2, line_segment.y2);
            result.emplace_back(line_segment);
        }
    }
    return result;
}

#ifdef VISUALIZATION
GrayScaledImage BuildEdgeMap(const std::vector<LineSegment>& line_segments, int32_t width, int32_t height) {
    GrayScaledImage result(width, height);
    for (const auto& line_segment : line_segments) {
        DrawLine(
                line_segment.y1,
                line_segment.x1,
                line_segment.y2,
                line_segment.x2,
                result
        );
    }
    return result;
}
#endif

Point Project(const Point& point, const Line& line) {
    double A = cos(line.angle);
    double B = sin(line.angle);
    double C = -line.r;
    double x = (B * (B * point.x - A * point.y) - A * C) / (A * A + B * B);
    double y = (A * (-B * point.x + A * point.y) - B * C) / (A * A + B * B);

    return Point{static_cast<int32_t>(y), static_cast<int32_t>(x)};
}

std::vector<Line> GetLines(std::vector<LineSegment>& line_segments, USE_PARAMS) {
    std::vector<Line> result;
    BoundedLineHeap heap(AMOUNT_OF_LINES_TO_PROCESS);
    std::sort(line_segments.begin(), line_segments.end(),
              [] (const LineSegment& l1, const LineSegment& l2) {
                  return l1.length > l2.length;
              });
    std::vector<bool> used(line_segments.size());
    ComparatorParam param;
    auto cmp = [&param] (const std::pair<Point, Point>& seg1,
                         const std::pair<Point, Point>& seg2) {
        if (param == ComparatorParam::ByX) {
            return seg1.first.x < seg2.first.x;
        }
        return seg1.first.y < seg2.first.y;
    };
    for (int32_t i = 0; i < line_segments.size(); ++i) {
        if (used[i]) {
            continue;
        }
        used[i] = true;
        Line line;
        line.weight = line_segments[i].length;
        line.angle = line_segments[i].angle - M_PI_2;
        if (line.angle < -M_PI) {
            line.angle += M2_PI;
        }
        line.r = ((line_segments[i].x1 + line_segments[i].x2) * cos(line.angle)
                  + (line_segments[i].y1 + line_segments[i].y2) * sin(line.angle)) / 2;
        line.segments.emplace_back(Point{line_segments[i].y1, line_segments[i].x1},
                                   Point{line_segments[i].y2, line_segments[i].x2});

        param = (fabs(cos(line.angle)) < fabs(sin(line.angle)))
                ? ComparatorParam::ByX : ComparatorParam::ByY;

        if (param == ComparatorParam::ByX) {
            if (line.segments.back().first.x > line.segments.back().second.x) {
                std::swap(line.segments.back().first, line.segments.back().second);
            }
        } else {
            if (line.segments.back().first.y > line.segments.back().second.y) {
                std::swap(line.segments.back().first, line.segments.back().second);
            }
        }

        for (int32_t j = i + 1; j < line_segments.size(); ++j) {
            double jangle = line_segments[j].angle - M_PI_2;
            if (jangle < -M_PI) {
                jangle += M2_PI;
            }
            double jr = ((line_segments[j].x1 + line_segments[j].x2) * cos(jangle)
                         + (line_segments[j].y1 + line_segments[j].y2) * sin(jangle)) / 2;
            if (IsAligned(line.angle, jangle,ANGLE_TOLERANCE)
                && fabs(line.r - jr) < LINE_R_TOLERANCE) {

                line.r = (line.weight * line.r + line_segments[j].length * jr)
                         / (line.weight + line_segments[j].length);
                double s = line.weight * sin(line.angle) + line_segments[j].length * sin(jangle)
                                                           / (line.weight + line_segments[j].length);
                double c = line.weight * cos(line.angle) + line_segments[j].length * cos(jangle)
                                                           / (line.weight + line_segments[j].length);
                line.angle = atan2(s, c);
                line.weight += line_segments[j].length;
                line.segments.emplace_back(Point{line_segments[j].y1, line_segments[j].x1},
                                           Point{line_segments[j].y2, line_segments[j].x2});

                if (param == ComparatorParam::ByX) {
                    if (line.segments.back().first.x > line.segments.back().second.x) {
                        std::swap(line.segments.back().first, line.segments.back().second);
                    }
                } else {
                    if (line.segments.back().first.y > line.segments.back().second.y) {
                        std::swap(line.segments.back().first, line.segments.back().second);
                    }
                }

                used[j] = true;
            }
        }
        for (int32_t k = 0; k < line.segments.size(); ++k) {
            line.segments[k] = std::make_pair(
                    Project(line.segments[k].first, line),
                    Project(line.segments[k].second, line)
            );
        }
        std::sort(line.segments.begin(), line.segments.end(), cmp);
        std::vector<std::pair<Point, Point>> collapsed_segments;
        std::pair<Point, Point> cur_collapsed_seg = line.segments[0];
        if (param == ComparatorParam::ByX) {
            for (int32_t k = 1; k < line.segments.size(); ++k) {
                if (cur_collapsed_seg.second.x >= line.segments[k].first.x) {
                    if (line.segments[k].second.x > cur_collapsed_seg.second.x) {
                        cur_collapsed_seg.second = line.segments[k].second;
                    }
                } else {
                    collapsed_segments.emplace_back(cur_collapsed_seg);
                    cur_collapsed_seg = line.segments[k];
                }
            }
        } else {
            for (int32_t k = 1; k < line.segments.size(); ++k) {
                if (cur_collapsed_seg.second.y >= line.segments[k].first.y) {
                    if (line.segments[k].second.y > cur_collapsed_seg.second.y) {
                        cur_collapsed_seg.second = line.segments[k].second;
                    }
                } else {
                    collapsed_segments.emplace_back(cur_collapsed_seg);
                    cur_collapsed_seg = line.segments[k];
                }
            }
        }
        collapsed_segments.emplace_back(cur_collapsed_seg);
        line.segments = std::move(collapsed_segments);
        line.weight = 0.0;
        for (auto& [p1, p2] : line.segments) {
            line.weight += Distance(p1.x, p1.y, p2.x, p2.y);
        }
        heap.Push(std::move(line));
    }
    while (!heap.Empty()) {
        result.emplace_back(heap.Pop());
    }
    std::reverse(result.begin(), result.end());
    return result;
}
