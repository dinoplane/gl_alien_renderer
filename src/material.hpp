#ifndef MATERIAL_H
#define MATERIAL_H

#include <glad/glad.h>
#include <glm/glm.hpp>

enum MaterialUniformFlags : std::uint32_t {
    None = 0 << 0,
    HasBaseColorTexture = 1 << 0,
};

struct Texture {
    GLuint texture;
};

struct Material {
    glm::vec4 baseColorFactor;
    float alphaCutoff;
    uint32_t flags;
};


#endif