cmake_minimum_required(VERSION 3.10)
project(_)

# 设置 CXX 版本
# 设置可执行文件输出路径
# 设置动态库文件输出路径
set(CMAKE_CXX_STANDARD 20)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/bin)
set(LIBRARY_OUTPUT_PATH ${PROJECT_SOURCE_DIR}/lib)

add_subdirectory(MultiThreadPrinter)
add_subdirectory(ImageFeatureExtractor)
add_subdirectory(SVGDSimulator)