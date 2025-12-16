find_package(OpenCV 4.7.0 REQUIRED)

if(OpenCV_FOUND)
    # After successfully finding the package
    message(STATUS "OpenCV 4.7.0 found successfully!")
    message(STATUS "OpenCV include directories: ${OpenCV_INCLUDE_DIRS}")
    message(STATUS "OpenCV libraries to link: ${OpenCV_LIBS}")
else()
    # This message is triggered by REQUIRED above, but good practice to show this explicitly:
    message(FATAL_ERROR "OpenCV 4.7.0 not found.")
endif()

add_executable(TestImageFeatureExtractor ImageFeatureExtractor_t.cpp)
