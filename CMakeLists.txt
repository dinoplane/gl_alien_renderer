cmake_minimum_required(VERSION 3.17)
project(flocking)

# Variables
set(VENDOR_DIR "${CMAKE_CURRENT_SOURCE_DIR}/vendor")
set(SRC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src")
set(SHADER_DIR "${CMAKE_CURRENT_SOURCE_DIR}/resources/shader")


set_property(GLOBAL PROPERTY USE_FOLDERS ON)
set (CMAKE_VERBOSE_MAKEFILE 0) # 1 should be used for debugging
set (CMAKE_SUPPRESS_REGENERATION TRUE) # Suppresses ZERO_CHECK
if(MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /W4")
else()
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wpedantic -std=c++11")
  if(NOT WIN32)
    set(GLAD_LIBRARIES dl)
  endif()
endif()

# OpenGL Flags
set(OpenGL_GL_PREFERENCE GLVND)

# GLFW - Installed via: sudo dnf install -y glfw-devel
option(GLFW_BUILD_DOCS OFF)
option(GLFW_BUILD_EXAMPLES OFF)
option(GLFW_BUILD_TESTS OFF)
option(GLFW_INSTALL OFF)
find_package(glfw3 REQUIRED)


#fmt
add_subdirectory(${VENDOR_DIR}/fmt EXCLUDE_FROM_ALL)


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
add_library(glm INTERFACE)
target_include_directories(glm INTERFACE ${VENDOR_DIR}/glm/glm)


# imgui
add_library(imgui STATIC)
target_include_directories(imgui PUBLIC ${VENDOR_DIR}/imgui)
target_sources(imgui PRIVATE
  ${VENDOR_DIR}/imgui/imconfig.h
  ${VENDOR_DIR}/imgui/imgui.h
  ${VENDOR_DIR}/imgui/imgui.cpp
  ${VENDOR_DIR}/imgui/imgui_demo.cpp
  ${VENDOR_DIR}/imgui/imgui_draw.cpp
  ${VENDOR_DIR}/imgui/imgui_widgets.cpp
  ${VENDOR_DIR}/imgui/imgui_tables.cpp
  ${VENDOR_DIR}/imgui/imgui_impl_glfw.cpp
  ${VENDOR_DIR}/imgui/imgui_impl_opengl3.cpp
  )


# stb
add_library (stb INTERFACE)
target_include_directories(stb INTERFACE ${VENDOR_DIR}/stb)

# Organize files
include_directories(
   ${VENDOR_DIR}/glm/
   ${VENDOR_DIR}/glad/include/
   ${VENDOR_DIR}/stb/
   ${SRC_DIR}/solids/
   ${SRC_DIR}/boids/
   ${SRC_DIR}/spatial_map/
   ${SRC_DIR}/
)


file(GLOB         VENDORS_SOURCES ${VENDOR_DIR}/glad/src/glad.c)
file(GLOB_RECURSE PROJECT_HEADERS ${SRC_DIR}/*.hpp
                                  ${SRC_DIR}/*.h)

file(GLOB_RECURSE PROJECT_SOURCES ${SRC_DIR}/*.cpp
                                  ${SRC_DIR}/*.cxx
                                  ${SRC_DIR}/*.cc
                                  ${SRC_DIR}/*.c)

file (GLOB_RECURSE PROJECT_SHADERS ${SHADER_DIR}/*.comp
                                   ${SHADER_DIR}/*.frag
                                   ${SHADER_DIR}/*.geom
                                   ${SHADER_DIR}/*.vert)

file (GLOB         PROJECT_CONFIGS CMakeLists.txt
                                 .gitignore
                                 .gitmodules)



# Project

source_group ("headers" FILES ${PROJECT_HEADERS})
source_group ("shaders" FILES ${PROJECT_SHADERS})
source_group ("sources" FILES ${PROJECT_SOURCES})
source_group ("vendors" FILES ${VENDORS_SOURCES})

add_definitions(-DGLFW_INCLUDE_NONE
                -DPROJECT_SOURCE_DIR=\"${PROJECT_SOURCE_DIR}\")
add_executable(${PROJECT_NAME} ${PROJECT_SOURCES} ${PROJECT_HEADERS}
                               ${PROJECT_SHADERS} ${PROJECT_CONFIGS}
                               ${VENDORS_SOURCES})

target_link_libraries(${PROJECT_NAME}
    PUBLIC
    glfw
    glad
    imgui
    ${CMAKE_DL_LIBS} # Needed for glad - https://stackoverflow.com/a/56842079/2394163
)

# set(CMAKE_BUILD_TYPE RelWithDebInfo)

set_target_properties(${PROJECT_NAME}
  PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
    COMPILE_OPTIONS -Wpedantic #-Werror
    CXX_STANDARD 20
    CXX_STANDARD_REQUIRED ON
    CXX_EXTENSIONS OFF
)