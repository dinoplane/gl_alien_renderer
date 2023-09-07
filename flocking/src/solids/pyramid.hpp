#ifndef PYRAMID_H
#define PYRAMID_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include<solid.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Pyramid : public Solid {
    public:
        virtual void init();
        virtual void calculateVertices();
};
#endif