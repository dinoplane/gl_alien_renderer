#ifndef MODEL_H
#define MODEL_H
#include <util.h>

#include<fastgltf/core.hpp>
#include <mesh.hpp>
#include <transform.hpp>
#include <volume.hpp>
#include <vector>
#include <glm/glm.hpp>

struct PrimitiveProperties {
    uint32_t materialIdx;
    uint32_t textureIdx;
};

struct NodeProperties {
    glm::mat4 modelFromMesh;
    // uint32_t primtiveIdx;
};

struct NodePrimProperties {
    uint32_t nodeIdx;
    uint32_t primIdx;
};

struct AtlasOffsets {
    GLsizei offsetX;
    GLsizei mipLevel; // Unused, but there for padding
    GLsizei textureWidth;
    GLsizei textureHeight;
};
struct GPUAtlasOffsets {
    GLfloat offsetX;
    GLfloat mipLevel; // Unused, but there for padding
    GLfloat textureWidth;
    GLfloat textureHeight;
};

struct Node {
    glm::mat4x4 nodeTransformMatrix;
    uint32_t meshIndex;

};

class Model {
    public:
    //  a model is a collection of meshes
    // each mesh can have its own material
    // each mesh can have its own texture
    // a model also contains a tree of nodes
    std::vector<Texture> textures;
    std::vector<Material> materials;

    // but in the end, a model could be a flat array of vertices and indices
    GLuint VBO, EBO;

    std::vector<IndirectDrawCommand> drawCmdBufferVec;
    GLuint drawCmdBuffer; // another compute pass to set the instance count in the buffer...?
    // prolly can try cpu side first in all honesty
    struct DrawCmdCount {
        uint32_t count;
    } drawCmdCount;
    GLuint drawCmdCountBuffer;

    std::vector<Primitive> primitives;
    std::vector<Mesh> meshes;
    std::vector<Node> nodes;

    std::vector<NodePrimProperties> nodePrimPropertiesBufferVec;
    GLuint nodePrimPropertiesBuffer;

    std::vector<NodeProperties> nodePropertiesBufferVec;
    GLuint nodePropertiesBuffer;

    std::vector<PrimitiveProperties> primitivePropertiesBufferVec;
    GLuint primitivePropertiesBuffer;

    GLuint materialPropertiesBuffer;
    GLuint textureArray;

    GLuint textureAtlas;
    std::vector<AtlasOffsets> atlasOffsetsBufferVec; // The first one will specify 0 and the atlas width
    std::vector<GPUAtlasOffsets> gpuAtlasOffsetsBufferVec;
    GLuint gpuAtlasOffsetsBuffer;

    // std::vector<glm::mat4>

/*

    // Struct for MultiDrawElements
    typedef  struct {
        unsigned int  count;         // index count
        unsigned int  instanceCount; // should always be 1
        unsigned int  firstIndex;
        int           baseVertex;
        unsigned int  baseInstance;
        // Optional user-defined data goes here - if nothing, stride is 0
    } DrawElementsIndirectCommand;

    void glMultiDrawElementsIndirect(
        GLenum mode,
        GLenum type,
        const void *indirect,
        GLsizei drawcount,
        GLsizei stride
    );


    // Goal: 1 draw call per model
    // Every model has its own materials/textures/primitives/nodes

    // flat array of GLuints of primitives in meshes

    In the renderer, it should look like this
    for ( const auto& [classname, entInstData] : scene.entityInstanceMap ){
        // do culling
        // bind textures
        // bind materials
        // draw from the indirect buffer


    }

 */


    Sphere boundingVolume;

    // void GenerateBoundingVolumeFromNodes();
    static Model CreateCube();
    static Model CreatePyramid();

};

#endif