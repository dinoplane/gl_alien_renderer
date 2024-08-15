#ifndef UTIL_H
#define UTIL_H

#include <glm/glm.hpp>

#define ARRAY_COUNT(x) (sizeof(x) / sizeof((x)[0]))

template <typename T, size_t N>
constexpr size_t array_size(T (&)[N]) {
    return N;
}


struct Vertex{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 texcoords;
};


#endif