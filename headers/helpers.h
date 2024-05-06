#pragma once

#include "image.h"

#include <cmath>
#include <cfloat>
#include <queue>

#define M2_PI 6.28318530718
#define RELATIVE_ERROR_FACTOR 100.0
#define ENOUGH_BIG 10'000'000
#define POINT_DEFAULT_RADIUS 5
#define HISTOGRAM_SIZE 256

class Histogram {
public:
    Histogram();

    void Inc(int32_t index);

    int32_t operator[](int32_t index) const;

    int32_t TotalCount() const;

private:
    int32_t total_count_ = 0;
    int32_t data_[HISTOGRAM_SIZE] = {0};
};

double CalculateHistogramMean(const Histogram& hist);

double CalculateHistogramVariance(const Histogram& hist);

enum ComparatorParam {
    ByX,
    ByY
};

struct Point {
    int32_t y = 0;
    int32_t x = 0;
};

struct LineSegment {
    int32_t x1 = 0;
    int32_t x2 = 0;
    int32_t y1 = 0;
    int32_t y2 = 0;
    int32_t length = 0;
    double angle = 0.0;
};

struct Line {
    double angle;
    double r;
    double weight;
    std::vector<std::pair<Point, Point>> segments;
};

class BoundedLineHeap {
public:
    BoundedLineHeap() = delete;

    BoundedLineHeap(int32_t size);

    void Push(Line line);

    Line Pop();

    int32_t Size() const;

    bool Empty() const;

private:
    struct LineComparator {
        bool operator()(const Line& line1, const Line& line2) const {
            return line1.weight > line2.weight;
        }
    };

    int32_t max_size_ = 0;
    std::priority_queue<Line, std::vector<Line>, LineComparator> heap_;
};

class Matrix {
public:
    Matrix() = delete;

    Matrix(int32_t size, int32_t value);

    int32_t* operator[](int32_t index);

    int32_t Size() const;

    Matrix(const Matrix& other) = delete;

    Matrix(Matrix&& other) = delete;

    Matrix operator=(const Matrix& other) = delete;

    Matrix operator=(Matrix&& other) = delete;

    ~Matrix();

private:
    int32_t size_ = 0;
    int32_t* data_ = nullptr;
};

Point Intersect(const Line& line1, const Line& line2);

static bool DoubleEqual(double a, double b) {
    double abs_diff,aa,bb,abs_max;
    if(a == b) {
        return true;
    }
    abs_diff = fabs(a - b);
    aa = fabs(a);
    bb = fabs(b);
    abs_max = aa > bb ? aa : bb;
    if(abs_max < DBL_MIN) abs_max = DBL_MIN;
    return (abs_diff / abs_max) <= (RELATIVE_ERROR_FACTOR * DBL_EPSILON);
}

static double Distance(double x1, double y1, double x2, double y2) {
    return sqrt((x1 - x2) * (x1 - x2) + (y2 - y1) * (y2 - y1));
}

static double Distance(Point p1, Point p2) {
    return sqrt((p1.x - p2.x) * (p1.x - p2.x) + (p1.y - p2.y) * (p1.y - p2.y));
}

struct Rectangle {
    double x1 = 0.0;
    double y1 = 0.0;
    double x2 = 0.0;
    double y2 = 0.0;
    double center_x = 0.0;
    double center_y = 0.0;
    double angle = 0.0;
    double width = 0.0;
    double vx = 0.0;
    double vy = 0.0;
    double probability = 0.0;
    double angle_tolerance = 0.0;
};

class RectangleIterator {
public:
    explicit RectangleIterator(const Rectangle& rectangle);

    Point GetPoint() const;

    void Next();

    bool IsEnd() const;

private:
    void Fill(int32_t y1, int32_t x1, int32_t y2, int32_t x2);

    int32_t min_x_;
    std::vector<int32_t> mins_;
    std::vector<int32_t> maxs_;
    int32_t index_;
    int32_t y_;
    int32_t x_;
};

int32_t CrossProduct(Point A, Point B, Point C);

bool IsConvex(Point p1, Point p2, Point p3, Point p4);

bool IsNotSelfIntersected(Point p1, Point p2, Point p3, Point p4);

double GetArea(const std::vector<Point>& points);

#ifdef VISUALIZATION
void DrawPoint(Point point, GrayScaledImage& img, int32_t color = WHITE, int32_t radius = POINT_DEFAULT_RADIUS);

void DrawLine(int32_t y1, int32_t x1, int32_t y2, int32_t x2, GrayScaledImage& img, int32_t color = WHITE);

void DrawPoint(Point point, RGBImage& img, RGB pixel = GREEN, int32_t radius = POINT_DEFAULT_RADIUS);

void DrawLine(int32_t y1, int32_t x1, int32_t y2, int32_t x2, RGBImage& img, RGB pixel = GREEN);

void DrawLine(const Line& line, GrayScaledImage& img);

GrayScaledImage DrawLines(const std::vector<Line>& lines, int32_t width, int32_t height);

void DrawCornerDetectorResult(const std::vector<Point>& corners, RGBImage& img);
#endif