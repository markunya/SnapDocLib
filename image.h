#pragma once

#include <cstdint>
#include <vector>

struct RGB {
    int32_t r, g, b;
};

template<typename T>
class Image {
public:
    Image() = default;

    Image(int32_t width, int32_t height) : width_(width), height_(height) {
        data_.resize(width * height);
    }

    int32_t Width() const {
        return width_;
    }

    int32_t Height() const {
        return height_;
    }

    void SetPixel(int32_t y, int32_t x, T pixel) {
        if (y < 0 || x < 0 || y >= height_ || x >= width_) {
            return;
        }
        data_[y * width_ + x] = pixel;
    }

    T GetPixel(int32_t y, int32_t x) const {
        if (x < 0) {
            x = 0;
        }
        if (x >= width_) {
            x = width_ - 1;
        }
        if (y < 0) {
            y = 0;
        }
        if (y >= height_) {
            y = height_ - 1;
        }
        return data_[y * width_ + x];
    }

    void Clear() {
        height_ = 0;
        width_ = 0;
        data_.clear();
    }

private:
    int32_t height_ = 0;
    int32_t width_ = 0;
    std::vector<T> data_;
};

using RGBImage = Image<RGB>;
using GrayScaledImage = Image<int32_t>;
using BoolImage = Image<bool>;
using DoubleImage = Image<double>;

inline int32_t ToGrayScale(const RGB &pixel) {
    return static_cast<int32_t >(0.299 * static_cast<double>(pixel.r) +
                                 0.587 * static_cast<double>(pixel.g) +
                                 0.114 * static_cast<double>(pixel.b));
}

inline GrayScaledImage RGBImageToGrayScaledImage(const RGBImage& img) {
    GrayScaledImage result(img.Width(), img.Height());
    for (int32_t i = 0; i < result.Height(); ++i) {
        for (int32_t j = 0; j < result.Width(); ++j) {
            result.SetPixel(i, j, ToGrayScale(img.GetPixel(i, j)));
        }
    }
    return result;
}

inline RGBImage GrayScaledImageToRGBImage(const GrayScaledImage& img) {
    RGBImage result(img.Width(), img.Height());
    for (int32_t i = 0; i < result.Height(); ++i) {
        for (int32_t j = 0; j < result.Width(); ++j) {
            result.SetPixel(i, j, {
                img.GetPixel(i, j),
                img.GetPixel(i, j),
                img.GetPixel(i, j),
            });
        }
    }
    return result;
}
