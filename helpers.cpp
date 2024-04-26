#include "helpers.h"

Histogram::Histogram() = default;

int32_t Histogram::operator[](int32_t index) const {
    return data_[index];
}

int32_t Histogram::TotalCount() const {
    return total_count_;
}

void Histogram::Inc(int32_t index) {
    ++data_[index];
    ++total_count_;
}


double CalculateHistogramMean(const Histogram& hist) {
    double sum = 0;
    int32_t total_count = 0;
    for (int32_t i = 0; i < HISTOGRAM_SIZE; ++i) {
        sum += i * hist[i];
        total_count += hist[i];
    }
    return total_count ? (sum / total_count) : 0;
}

double CalculateHistogramVariance(const Histogram& hist) {
    double mean = CalculateHistogramMean(hist);
    double variance = 0;
    int32_t total_count = hist.TotalCount();

    for (int i = 0; i < 256; ++i) {
        double diff = i - mean;
        variance += diff * diff * hist[i];
    }

    return total_count ? (variance / total_count) : 0;
}

Point Intersect(const Line& line1, const Line& line2) {
    double A1 = cos(line1.angle);
    double B1 = sin(line1.angle);
    double C1 = -line1.r;

    double A2 = cos(line2.angle);
    double B2 = sin(line2.angle);
    double C2 = -line2.r;

    double denominator = A1*B2 - A2*B1;
    if (fabs(denominator) < 1e-7) {
        return Point{ENOUGH_BIG, ENOUGH_BIG};
    }

    double x = (B1*C2 - B2*C1) / denominator;
    double y = (C1*A2 - C2*A1) / denominator;
    if (fabs(x) < ENOUGH_BIG && fabs(y) < ENOUGH_BIG) {
        return Point{static_cast<int32_t>(std::round(y)), static_cast<int32_t>(std::round(x))};
    }
    return Point{ENOUGH_BIG, ENOUGH_BIG};
}

RectangleIterator::RectangleIterator(const Rectangle &rectangle) {
    int32_t corner1_x, corner2_x, corner3_x, corner4_x;
    int32_t corner1_y, corner2_y, corner3_y, corner4_y;
    double bayes_x = (rectangle.width / 2) * cos(rectangle.angle + (M_PI / 2));
    double bayes_y = (rectangle.width / 2) * sin(rectangle.angle + (M_PI / 2));
    corner1_x = static_cast<int32_t>(rectangle.x1 + bayes_x);
    corner2_x = static_cast<int32_t>(rectangle.x2 + bayes_x);
    corner3_x = static_cast<int32_t>(rectangle.x2 - bayes_x);
    corner4_x = static_cast<int32_t>(rectangle.x1 - bayes_x);
    corner1_y = static_cast<int32_t>(rectangle.y1 + bayes_y);
    corner2_y = static_cast<int32_t>(rectangle.y2 + bayes_y);
    corner3_y = static_cast<int32_t>(rectangle.y2 - bayes_y);
    corner4_y = static_cast<int32_t>(rectangle.y1 - bayes_y);
    int32_t min_x = std::min(corner1_x, std::min(corner2_x, std::min(corner3_x, corner4_x)));
    int32_t max_x = std::max(corner1_x, std::max(corner2_x, std::max(corner3_x, corner4_x)));
    mins_.resize(max_x - min_x + 1, std::numeric_limits<int32_t>::max());
    maxs_.resize(max_x - min_x + 1, std::numeric_limits<int32_t>::min());
    min_x_ = min_x;
    Fill(corner1_y, corner1_x, corner2_y, corner2_x);
    Fill(corner2_y, corner2_x, corner3_y, corner3_x);
    Fill(corner3_y, corner3_x, corner4_y, corner4_x);
    Fill(corner4_y, corner4_x, corner1_y, corner1_x);
    index_ = 0;
    x_ = min_x_;
    y_ = maxs_[0];
}

Point RectangleIterator::GetPoint() const{
    return Point{y_, x_};
}

void RectangleIterator::Next() {
    if (index_ == mins_.size()) {
        return;
    }
    if (y_ > mins_[index_]) {
        --y_;
        return;
    }
    ++index_;
    if (index_ == maxs_.size()) {
        return;
    }
    ++x_;
    y_ = maxs_[index_];
}

bool RectangleIterator::IsEnd() const {
    return index_ == mins_.size();
}

void RectangleIterator::Fill(int32_t y1, int32_t x1, int32_t y2, int32_t x2) {
    int32_t dx = std::abs(x2 - x1);
    int32_t dy = std::abs(y2 - y1);
    int32_t sx = x1 < x2 ? 1 : -1;
    int32_t sy = y1 < y2 ? 1 : -1;
    int error = dx - dy;
    mins_[x2 - min_x_] = std::min(y2, mins_[x2 - min_x_]);
    maxs_[x2 - min_x_] = std::max(y2, maxs_[x2 - min_x_]);
    while(x1 != x2 || y1 != y2) {
        mins_[x1 - min_x_] = std::min(y1, mins_[x1 - min_x_]);
        maxs_[x1 - min_x_] = std::max(y1, maxs_[x1 - min_x_]);
        int error2 = error * 2;
        if(error2 > -dy) {
            error -= dy;
            x1 += sx;
        }
        if(error2 < dx) {
            error += dx;
            y1 += sy;
        }
    }
}

Matrix::Matrix(int32_t size, int32_t value) : size_(size) {
    data_ = new int32_t[size * size];
    for (int32_t i = 0; i < size * size; ++i) {
        data_[i] = value;
    }
}

int32_t *Matrix::operator[](int32_t index) {
    return data_ + index * size_;
}

Matrix::~Matrix() {
    delete[] data_;
}

int32_t Matrix::Size() const {
    return size_;
}

int32_t CrossProduct(Point A, Point B, Point C) {
    return (B.x - A.x) * (C.y - B.y) - (B.y - A.y) * (C.x - B.x);
}

bool IsConvex(Point p1, Point p2, Point p3, Point p4) {
    int32_t cp1 = CrossProduct(p1, p2, p3);
    int32_t cp2 = CrossProduct(p2, p3, p4);
    int32_t cp3 = CrossProduct(p3, p4, p1);
    int32_t cp4 = CrossProduct(p4, p1, p2);

    if ((cp1 >= 0 && cp2 >= 0 && cp3 >= 0 && cp4 >= 0) ||
        (cp1 <= 0 && cp2 <= 0 && cp3 <= 0 && cp4 <= 0)) {
        return true;
    }
    return false;
}

bool IsNotSelfIntersected(Point p1, Point p2, Point p3, Point p4) {
    int32_t cp1 = CrossProduct(p2, p1, p4);
    int32_t cp2 = CrossProduct(p2, p4, p3);
    if ((cp1 >= 0 && cp2 >= 0) || (cp1 <= 0 && cp2 <= 0)) {
        return true;
    }
    return false;
}

double GetArea(const std::vector<Point>& points) {
    double area = 0;
    for (int32_t i = 0; i < points.size(); ++i) {
        area += (points[(i + 1) % points.size()].x - points[i].x)
                * (static_cast<double>((points[i].y + points[(i + 1) % points.size()].y)) / 2);
    }
    return fabs(area);
}

BoundedLineHeap::BoundedLineHeap(int32_t size) : max_size_(size) {
}

void BoundedLineHeap::Push(Line line) {
    heap_.emplace(std::move(line));
    if (heap_.size() > max_size_) {
        heap_.pop();
    }
}

Line BoundedLineHeap::Pop() {
    Line l = heap_.top();
    heap_.pop();
    return l;
}

bool BoundedLineHeap::Empty() const {
    return heap_.empty();
}

#ifdef VISUALIZATION
void DrawLine(int32_t y1, int32_t x1, int32_t y2, int32_t x2, GrayScaledImage& img, int32_t color) {
    int32_t dx = std::abs(x2 - x1);
    int32_t dy = std::abs(y2 - y1);
    int32_t sx = x1 < x2 ? 1 : -1;
    int32_t sy = y1 < y2 ? 1 : -1;
    int error = dx - dy;
    img.SetPixel(y2, x2, color);
    while(x1 != x2 || y1 != y2) {
        img.SetPixel(y1, x1, color);
        int error2 = error * 2;
        if(error2 > -dy) {
            error -= dy;
            x1 += sx;
        }
        if(error2 < dx) {
            error += dx;
            y1 += sy;
        }
    }
}

void DrawPoint(Point point, GrayScaledImage& img, int32_t color, int32_t radius) {
    for (int32_t i = -radius; i <= radius; ++i) {
        for (int32_t j = -radius; j <= radius; ++j) {
            img.SetPixel(point.y + i, point.x + j, color);
        }
    }
}

void DrawPoint(Point point, RGBImage& img, RGB color, int32_t radius) {
    for (int32_t i = -radius; i <= radius; ++i) {
        for (int32_t j = -radius; j <= radius; ++j) {
            img.SetPixel(point.y + i, point.x + j, color);
        }
    }
}

void DrawLine(int32_t y1, int32_t x1, int32_t y2, int32_t x2, RGBImage& img, RGB color) {
    int32_t dx = std::abs(x2 - x1);
    int32_t dy = std::abs(y2 - y1);
    int32_t sx = x1 < x2 ? 1 : -1;
    int32_t sy = y1 < y2 ? 1 : -1;
    int error = dx - dy;
    img.SetPixel(y2, x2, color);
    while(x1 != x2 || y1 != y2) {
        img.SetPixel(y1, x1, color);
        int error2 = error * 2;
        if(error2 > -dy) {
            error -= dy;
            x1 += sx;
        }
        if(error2 < dx) {
            error += dx;
            y1 += sy;
        }
    }
}

void DrawLine(const Line& line, GrayScaledImage& img) {
    //r = x * cos(a) + y * sin(a)
    double s = sin(line.angle);
    double c = cos(line.angle);
    double x1, y1, x2, y2;
    if (fabs(c) < fabs(s)) {
        // y = r / sin(a) - x * ctg(a)
        x1 = 0.0;
        x2 = img.Width();
        y1 = line.r / s - (x1 * c) / s;
        y2 = line.r / s - (x2 * c) / s;
    } else {
        // x = r / cos(a) - y * tg(a)
        y1 = 0.0;
        y2 = img.Height();
        x1 = line.r / c - (y1 * s) / c;
        x2 = line.r / c - (y2 * s) / c;
    }
    DrawLine(
            static_cast<int32_t>(std::round(y1)),
            static_cast<int32_t>(std::round(x1)),
            static_cast<int32_t>(std::round(y2)),
            static_cast<int32_t>(std::round(x2)), img, GRAY);
    for (const auto& [p1, p2] : line.segments) {
        DrawLine(p1.y, p1.x, p2.y, p2.x, img, WHITE);
    }
}

GrayScaledImage DrawLines(const std::vector<Line>& lines, int32_t width, int32_t height) {
    GrayScaledImage result(width, height);
    for (const auto& line : lines) {
        //r = x * cos(a) + y * sin(a)
        double s = sin(line.angle);
        double c = cos(line.angle);
        double x1, y1, x2, y2;
        if (fabs(c) < fabs(s)) {
            // y = r / sin(a) - x * ctg(a)
            x1 = 0.0;
            x2 = width;
            y1 = line.r / s - (x1 * c) / s;
            y2 = line.r / s - (x2 * c) / s;
        } else {
            // x = r / cos(a) - y * tg(a)
            y1 = 0.0;
            y2 = height;
            x1 = line.r / c - (y1 * s) / c;
            x2 = line.r / c - (y2 * s) / c;
        }
        DrawLine(
                static_cast<int32_t>(std::round(y1)),
                static_cast<int32_t>(std::round(x1)),
                static_cast<int32_t>(std::round(y2)),
                static_cast<int32_t>(std::round(x2)), result, GRAY);
    }
    for (const auto& line : lines) {
        for (const auto& [p1, p2] : line.segments) {
            DrawLine(p1.y, p1.x, p2.y, p2.x, result, WHITE);
        }
    }
    return result;
}

void DrawCornerDetectorResult(const std::vector<Point>& corners, RGBImage& img) {
    for (int32_t i = 0; i < corners.size(); ++i) {
        DrawLine(
                corners[i].y,
                corners[i].x,
                corners[(i + 1) % corners.size()].y,
                corners[(i + 1) % corners.size()].x,
                img
        );
    }
    for (int32_t i = 0; i < corners.size(); ++i) {
        DrawPoint(corners[i], img);
    }
}
#endif
