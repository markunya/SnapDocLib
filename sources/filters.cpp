#include "../headers/filters.h"

#include <cmath>
#include <algorithm>

GrayScaledImage GaussianBlur(const GrayScaledImage& img) {
    int32_t normalize_k = 159;
    std::vector<std::vector<int32_t>> matrix = {
            {2, 4, 5, 4, 2},
            {4, 9, 12, 9, 4},
            {5, 12, 15, 12, 5},
            {4, 9, 12, 9, 4},
            {2, 4, 5, 4, 2},
    };
    int32_t bayes = static_cast<int32_t>(matrix.size()) / 2;
    GrayScaledImage result(img.Width(), img.Height());
    for (int32_t i = 0; i < img.Height(); ++i) {
        for (int32_t j = 0; j < img.Width(); ++j) {
            int pixel = 0;
            for (int32_t k = -bayes; k <= bayes; ++k) {
                for (int32_t l = -bayes; l <= bayes; ++l) {
                    if (i + k >= 0 && i + k < static_cast<int32_t>(img.Height())
                    && j + l >= 0 && j + l < static_cast<int32_t>(img.Width())) {
                        pixel += matrix[k + bayes][l + bayes] * img.GetPixel(i + k, j + l);
                    } else {
                        pixel += matrix[k + bayes][l + bayes] * img.GetPixel(i, j);
                    }
                }
            }
            pixel /= normalize_k;
            result.SetPixel(i, j, pixel);
        }
    }
    return result;
}

GrayScaledImage Reduce(const GrayScaledImage& img, int32_t& reduce_coefficient, USE_PARAMS) {
    if (img.Width() * img.Height() <= MAX_SIZE) {
        return img;
    }
    auto width = static_cast<double>(img.Width());
    auto height = static_cast<double>(img.Height());
    reduce_coefficient = std::ceil(sqrt(width * height / MAX_SIZE));
    GrayScaledImage result(img.Width() / reduce_coefficient, img.Height() / reduce_coefficient);
    for (int32_t i = 0; i < height; ++i) {
        for (int32_t j = 0; j < width; ++j) {
            int32_t k = i / reduce_coefficient;
            int32_t l = j / reduce_coefficient;
            result.SetPixel(k, l, result.GetPixel(k, l) + img.GetPixel(i, j));
        }
    }
    for (int32_t i = 0; i < result.Height(); ++i) {
        for (int32_t j = 0; j < result.Width(); ++j) {
            result.SetPixel(i, j, result.GetPixel(i, j) / (reduce_coefficient * reduce_coefficient));
        }
    }
    return result;
}
