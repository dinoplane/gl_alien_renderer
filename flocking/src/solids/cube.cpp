#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include<cube.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

void Cube::init(){
    calculateVertices();
    initializeBuffers();
}

void Cube::calculateVertices(){
    std::cout << "Hello" << std::endl;
    vertices = {
        // positions         // colors              // texcoords

         0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,     1.0, 0.0,       // top right
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,     1.0, 1.0,       // bottom right
        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 0.0f,     0.0, 1.0,       // bottom left
        -0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 0.0f,     0.0, 0.0,       // top left
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 1.0f,     0.0, 0.0,       // top right
         0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 1.0f,     0.0, 1.0,       // bottom right
        -0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 1.0f,     1.0, 1.0,       // bottom left
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,     1.0, 0.0,       // top left

         0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,     0.0, 0.0,       // top right
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,     0.0, 1.0,       // bottom right
        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 0.0f,     1.0, 1.0,       // bottom left
        -0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 0.0f,     1.0, 0.0,       // top left
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 1.0f,     1.0, 0.0,       // top right
         0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 1.0f,     1.0, 1.0,       // bottom right
        -0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 1.0f,     0.0, 1.0,       // bottom left
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,     0.0, 0.0,       // top left

         0.5f,  0.5f,  0.5f,  1.0f, 1.0f, 0.0f,     1.0, 1.0,       // top right
         0.5f, -0.5f,  0.5f,  1.0f, 0.0f, 0.0f,     1.0, 0.0,       // bottom right
        -0.5f, -0.5f,  0.5f,  0.0f, 0.0f, 0.0f,     0.0, 0.0,       // bottom left
        -0.5f,  0.5f,  0.5f,  0.0f, 1.0f, 0.0f,     0.0, 1.0,       // top left
         0.5f,  0.5f, -0.5f,  1.0f, 1.0f, 1.0f,     1.0, 0.0,       // top right
         0.5f, -0.5f, -0.5f,  1.0f, 0.0f, 1.0f,     1.0, 1.0,       // bottom right
        -0.5f, -0.5f, -0.5f,  0.0f, 0.0f, 1.0f,     0.0, 1.0,       // bottom left
        -0.5f,  0.5f, -0.5f,  0.0f, 1.0f, 1.0f,     0.0, 0.0,       // top left
    };



    indices = {  // note that we start from 0!
        0, 1, 3,  // front tr br tl
        1, 2, 3,  // front br bl tl

        12, 13, 8,  // right tr br tl
        13, 9, 8,   // right br bl tl

        20, 16, 23,  // top tr br tl
        16, 19, 23,  // top br bl tl

        7, 6, 4,  // back tr br tl
        6, 5, 4,  // back br bl tl

        11, 10, 15,  // left tr br tl
        10, 14, 15,  // left br bl tl

        17, 21, 18,  // bot tr br tl
        21, 22, 18   // bot br bl tl
    };

}