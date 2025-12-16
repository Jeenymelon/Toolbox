# 1. Find the packages
find_package(OpenGL REQUIRED)
find_package(glfw3 CONFIG QUIET)

if(NOT glfw3_FOUND)
    find_package(PkgConfig REQUIRED)
    pkg_check_modules(GLFW REQUIRED glfw3)
endif()

# 2. Define your executable (You already have this line)
#add_executable(TestSVGDSimulator P2dSVGDSimulator_t.cpp)
add_executable(TestSVGDSimulator P3dSVGDSimulator_t.cpp)

# 3. LINK LIBRARIES (This is the missing part!)
# You must link OpenGL::GL and OpenGL::GLU to 'TestSVGDSimulator'
if(glfw3_FOUND)
    target_link_libraries(TestSVGDSimulator PRIVATE glfw OpenGL::GL OpenGL::GLU)
else()
    target_link_libraries(TestSVGDSimulator PRIVATE ${GLFW_LIBRARIES} OpenGL::GL OpenGL::GLU)
    target_include_directories(TestSVGDSimulator PRIVATE ${GLFW_INCLUDE_DIRS})
endif()