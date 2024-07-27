#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


#include <glm/gtx/string_cast.hpp>


#include<plane.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

void Plane::init(){
    calculateVertices();
    initializeBuffers();
}

void Plane::calculateVertices(){
    vertices.clear();
    indices.clear();
    // std::cout << "Hello" << std::endl;

    //  point of each (g+1)^2 vertices
    float vx = -0.5f;
    float vy = -0.5f;

    for (int j = 0; j <= grid_dims.y; j++){
        vx = -0.5f;
        for (int i = 0; i <= grid_dims.x; i++){
            float a[8] = {
                vx,  0.0f,  vy,    1.0f, 1.0f, 0.0f,   vx + 0.5f, vy + 0.5f
            };

            vertices.insert(vertices.end(), std::begin(a), std::end(a));     // C++11

            vx += unit_dims.x;
        }
        vy += unit_dims.y;
    }

    for (int j = 0; j < grid_dims.y; j++){
        for (int i = 0; i < grid_dims.x; i++){
            // (i, j)     (i+1, j) (i+1, j+1)
            // (i+1, j+1) (i, j+1) (i, j)

            int tl_i = j*(grid_dims.x + 1) + i;

            int b[6] = {
                tl_i,                   tl_i + 1,           tl_i + 1 + (grid_dims.x + 1),
                tl_i + 1 + (grid_dims.x + 1), tl_i + (grid_dims.x + 1), tl_i
            };

            indices.insert(indices.end(), std::begin(b), std::end(b));     // C++11
        }
    }

    // printVertices();
}


void Plane::setTiles(glm::ivec2 gdims)
{
    grid_dims = gdims;
    unit_dims = glm::vec2(1.0f) / glm::vec2(grid_dims);
    // std::cout << glm::to_string(grid_dims) << std::endl;
    // std::cout << glm::to_string(unit_dims) << std::endl;
};