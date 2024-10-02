#ifndef MATERIAL_H
#define MATERIAL_H

#include <glad/glad.h>
#include <glm/glm.hpp>

enum MaterialUniformFlags : uint {
    None = 0 << 0,
    HasBaseColorTexture = 1 << 0,
};

struct Texture {
    GLuint texture;
};


struct alignas(16) Material {
    glm::vec4 baseColorFactor;
    float alphaCutoff;
    uint flags;
};

// struct Material{
//     MaterialUniforms uniforms;
//     std::vector<GLuint> 
// }       

#endif