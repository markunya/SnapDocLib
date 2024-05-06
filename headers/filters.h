#pragma once

#include "image.h"
#include "cd_params.h"

#include <vector>
#include <cstdint>

GrayScaledImage GaussianBlur(const GrayScaledImage& img);
GrayScaledImage Reduce(const GrayScaledImage& img, int32_t& reduce_coefficient, USE_PARAMS);

