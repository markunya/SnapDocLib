#pragma once

#include <cmath>

#define INIT_PARAMS CornerDetectorParams params;
#define USE_PARAMS const CornerDetectorParams& params
#define PARAMS params
#define IMG_PATH params.path
#define MAX_SIZE params.max_size
#define AIRBAG params.airbag
#define LINE_R_TOLERANCE params.line_r_tolerance
#define AMOUNT_OF_LINES_TO_PROCESS params.amount_of_lines_to_process
#define OPTIONAL_SEGMENT_SIZE params.optional_segment_size
#define OPTIONAL_SEGMENT_PERIOD params.optional_segment_period
#define PD_K params.pd_k
#define AP_K params.ap_k
#define SP_K params.sp_k
#define MIN_REG_SIZE_K params.lsd_params.min_reg_size_k
#define QUANT params.lsd_params.quant
#define ANGLE_TOLERANCE params.lsd_params.angle_tolerance
#define PROBABILITY params.lsd_params.probability
#define REQ_DENSITY params.lsd_params.req_density
#define GRAD_MAGNITUDE_THRESHOLD params.lsd_params.gradient_magnitude_threshold
#define LOG_EPS params.lsd_params.log_eps

struct LSDParams {
    double min_reg_size_k = 5.0;
    double quant = 2.0;
    double angle_tolerance = M_PI / 8.0;
    double probability = angle_tolerance / M_PI;
    double req_density = 0.7;
    double gradient_magnitude_threshold = quant / sin(angle_tolerance);
    double log_eps = 0.0;
};

struct CornerDetectorParams {
    const char* path = nullptr;
    int32_t max_size = 40'000;
    double line_r_tolerance = 8.0;
    int32_t amount_of_lines_to_process = 25;
    int32_t airbag  = 30;
    int32_t optional_segment_size = 10;
    int32_t optional_segment_period = 5;
    double pd_k = 50.0;
    double ap_k = 25.0;
    double sp_k = 25.0;
    LSDParams lsd_params;
};

