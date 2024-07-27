#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include<pyramid.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

void Pyramid::init(){
    calculateVertices();
    initializeBuffers();
}

void Pyramid::calculateVertices(){
    // std::cout << "Hello" << std::endl;
    vertices = {
        // positions         // colors              // texcoords
         0.5f, -0.5f,  0.5f,    1.0f, 1.0f, 0.0f,   1.0f, 1.0f,
        -0.5f, -0.5f,  0.5f,    1.0f, 1.0f, 0.0f,   0.0f, 1.0f,
        -0.5f, -0.5f, -0.5f,    1.0f, 1.0f, 0.0f,   1.0f, 1.0f,
         0.5f, -0.5f, -0.5f,    1.0f, 1.0f, 0.0f,   0.0f, 1.0f,
         0.0f,  0.5f,  0.0f,    1.0f, 1.0f, 0.0f,   0.5f, 0.0f,

         0.5f, -0.5f,  0.5f,    1.0f, 1.0f, 0.0f,   0.0f, 1.0f,
        -0.5f, -0.5f,  0.5f,    1.0f, 1.0f, 0.0f,   1.0f, 1.0f,
        -0.5f, -0.5f, -0.5f,    1.0f, 1.0f, 0.0f,   0.0f, 1.0f,
         0.5f, -0.5f, -0.5f,    1.0f, 1.0f, 0.0f,   1.0f, 1.0f,

         0.5f, -0.5f,  0.5f,    1.0f, 1.0f, 0.0f,   1.0f, 0.0f,
        -0.5f, -0.5f,  0.5f,    1.0f, 1.0f, 0.0f,   0.0f, 0.0f,
        -0.5f, -0.5f, -0.5f,    1.0f, 1.0f, 0.0f,   0.0f, 1.0f,
         0.5f, -0.5f, -0.5f,    1.0f, 1.0f, 0.0f,   1.0f, 1.0f
    };



    indices = {  // note that we start from 0!
        0, 1, 4,  // front tr br tl

        5, 8, 4,  // right tr br tl

        2, 3, 4,  // back tr br tl

        6, 7, 4,  // left tr br tl

        9, 12, 10,  // bot tr br tl
        12, 11, 10   // bot br bl tl
    };

}