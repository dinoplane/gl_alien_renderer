#ifndef MODEL_H
#define MODEL_H
#include <util.h>

#include<fastgltf/core.hpp>
#include <mesh.hpp>
#include <transform.hpp>
#include <volume.hpp>
#include <vector>
#include <glm/glm.hpp>

struct Node {
    glm::mat4x4 nodeTransformMatrix;
    uint meshIndex;

};

class Model {
    public: 
    //  a model is a collection of meshes
    // each mesh can have its own material
    // each mesh can have its own texture
    // a model also contains a tree of nodes
    std::vector<Texture> textures;
    std::vector<Material> materials;
    std::vector<Mesh> meshes;
    std::vector<Node> nodes;
    
    Sphere boundingVolume;

    // void GenerateBoundingVolumeFromNodes();
    static Model CreateCube();
    static Model CreatePyramid();
    
};

#endif