#ifndef PLANE_H
#define PLANE_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include<solid.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Plane : public Solid {
    public:
        glm::ivec2 grid_dims = glm::ivec2(1);
        glm::vec2 unit_dims = glm::vec2(1.0f);

        virtual void init();
        virtual void calculateVertices();

        void setTiles(glm::ivec2 gdims){
            grid_dims = gdims;
            unit_dims = glm::vec2(1.0f) / glm::vec2(grid_dims);
            std::cout << glm::to_string(grid_dims) << std::endl;
            std::cout << glm::to_string(unit_dims) << std::endl;

        };
};
#endif