cmake_minimum_required(VERSION 3.28)
project(SnapDocLib)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake-modules)

set(OpenCV_DIR "/home/xubuntu/OpenCV-android-sdk/sdk/native/jni")

add_library(SnapDocLib SHARED
        snap_doc_lib.cpp
        candidate_discoverer.cpp
        filters.cpp
        helpers.cpp
        lsd.cpp
        snap_doc_lib.cpp)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})
target_link_libraries(SnapDocLib ${OpenCV_LIBS})
