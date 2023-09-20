cmake_minimum_required(VERSION 3.17)
project(flocking)

# Variables
set(VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/vendor")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/${PROJECT_NAME}/src")


# OpenGL Flags
set(OpenGL_GL_PREFERENCE GLVND)

# GLFW - Installed via: sudo dnf install -y glfw-devel
option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_INSTALL OFF)
find_package(glfw3 REQUIRED)


# Glad - Installed via: git submodule add git@github.com:Dav1dde/glad.git videogame/vendor/glad
add_library(glad
  STATIC
    ${VENDOR_DIR}/glad/src/glad.c
)

target_include_directories(glad
  PUBLIC
    ${VENDOR_DIR}/glad/include
)

# glm
# add_library(glm
#   STATIC
#     ${VENDOR_DIR}/glad/src/glad.c
# )

include_directories(
   ${VENDOR_DIR}/glm/
   ${VENDOR_DIR}/glad/include/
   ${VENDOR_DIR}/stb/
   ${SRC_DIR}/solids/
   ${SRC_DIR}/boids/
   ${SRC_DIR}/spatial_map/
   ${SRC_DIR}/
)


# file(GLOB VENDORS_SOURCES Glitter/Vendor/glad/src/glad.c)
file(GLOB PROJECT_HEADERS ${SRC_DIR}/*.hpp
                          ${SRC_DIR}/solids/*.hpp
                          ${SRC_DIR}/boids/*.hpp
                          ${SRC_DIR}/spatial_map/*.hpp

                          )

file(GLOB PROJECT_SOURCES ${SRC_DIR}/*.cpp
                          ${SRC_DIR}/solids/*.cpp
                          ${SRC_DIR}/boids/*.cpp
                          ${SRC_DIR}/spatial_map/*.cpp

                          )
# file(GLOB PROJECT_SHADERS Glitter/Shaders/*.comp
#                           Glitter/Shaders/*.frag
#                           Glitter/Shaders/*.geom
#                           Glitter/Shaders/*.vert)
# file(GLOB PROJECT_CONFIGS CMakeLists.txt
#                           Readme.md
#                          .gitattributes
#                          .gitignore
#                          .gitmodules)

source_group("Headers" FILES ${PROJECT_HEADERS})
# source_group("Shaders" FILES ${PROJECT_SHADERS})
source_group("Sources" FILES ${PROJECT_SOURCES})
# source_group("Vendors" FILES ${VENDORS_SOURCES})

add_definitions(-DGLFW_INCLUDE_NONE
                -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS})
                              #  ${PROJECT_SHADERS} ${PROJECT_CONFIGS}
                              #  ${VENDORS_SOURCES})

# Project
# add_executable(${PROJECT_NAME})

# target_include_directories(${PROJECT_NAME}
#   PRIVATE
#     ${VENDOR_DIR}/glad/include
#     ${SRC_DIR}/solids
#     ${SRC_DIR}/

# )

# target_sources(${PROJECT_NAME}
#   PRIVATE
#     ${SRC_DIR}/main.cpp
#     ${SRC_DIR}/shader_s.h
#     ${SRC_DIR}/solids/*.cpp
# )

target_link_libraries(${PROJECT_NAME}
  PRIVATE
    glfw
    glad
    ${CMAKE_DL_LIBS} # Needed for glad - https://stackoverflow.com/a/56842079/2394163
)

# set(CMAKE_BUILD_TYPE RelWithDebInfo)

set_target_properties(${PROJECT_NAME}
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/${PROJECT_NAME}
    COMPILE_OPTIONS -Werror
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)