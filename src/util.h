#ifndef UTIL_H
#define UTIL_H

#include <glm/glm.hpp>

#define ARRAY_COUNT(x) (sizeof(x) / sizeof((x)[0]))
struct Vertex{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 texcoords;
};


#endif