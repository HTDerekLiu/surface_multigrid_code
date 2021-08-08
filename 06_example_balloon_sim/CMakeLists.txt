cmake_minimum_required(VERSION 3.1)
SET(CMAKE_BUILD_TYPE release)
project(main)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

# libigl
option(LIBIGL_WITH_OPENGL            "Use OpenGL"         ON)
option(LIBIGL_WITH_OPENGL_GLFW       "Use GLFW"           ON)
option(LIBIGL_WITH_OPENGL_GLFW_IMGUI "Use ImGui"          ON)
option(LIBIGL_WITH_VIEWER            "Use OpenGL viewer"  ON)

find_package(LIBIGL REQUIRED)

# Add your project files
include_directories("../utils/" "./simulation_utils")

get_property(dirs DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR} PROPERTY INCLUDE_DIRECTORIES)
foreach(dir ${dirs})
  message(STATUS "dir='${dir}'")
endforeach()

file(GLOB SRCFILES *.cpp ../utils/*.cpp ./simulation_utils/*.cpp ./simulation_utils/SecondFundamentalForm/*.cpp ./simulation_utils/MaterialModel/*.cpp)
add_executable(${PROJECT_NAME}_bin ${SRCFILES})

target_link_libraries(${PROJECT_NAME}_bin igl::core igl::opengl igl::opengl_glfw igl::opengl_glfw_imgui)



