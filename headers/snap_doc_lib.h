#pragma once
#include <cinttypes>

extern "C" {

int32_t DetectCorners(const char* path, int32_t* corners_ptr);
int32_t Normalize(const char* in_path, const char* out_path, int32_t* corners_ptr);

}
