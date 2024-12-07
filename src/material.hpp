#ifndef MATERIAL_H
#define MATERIAL_H

#include <glad/glad.h>
#include <glm/glm.hpp>

enum MaterialUniformFlags : uint32_t {
    None = 0 << 0,
    HasBaseColorTexture = 1 << 0,
};

struct Texture {
    GLuint texture;
    GLsizei textureWidth;
    GLsizei textureHeight;
    GLsizei mipLevels;
};


struct alignas(16) Material {
    glm::vec4 baseColorFactor;
    float alphaCutoff;
    uint32_t flags;
};

// struct Material{
//     MaterialUniforms uniforms;
//     std::vector<GLuint> 
// }       

#endif