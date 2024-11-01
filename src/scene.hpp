#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <unordered_map>
// #include <memory>
#include <mesh.hpp>
#include <model.hpp>
#include <gpu_structs.hpp>
class Entity;
class Shader; // Change this to only use move semantics
class Camera;
// class Light;
// class Material;
// class Entity;

struct EntityInstanceData
{
    Model instModel; // essentially the bounding volumes and the mesh only need one of...
    uint instCount;
    std::vector<glm::mat4> modelToWorldMat;
    GLuint instModelMatrixBuffer;
    std::vector<int> isInstMeshRendered;
    GLuint instMeshRenderedBuffer;

    // struct VisibleInstIndicesBlock {
    //     uint visibleInstCount;
    // //     std::vector<uint> visibleInstIndices;
    // } visibleInstIndicesBlock;

    std::vector<uint> visibleInstIndicesBlock;
    GLuint visibleInstIndicesSSBO;

    void GenerateInstanceBuffers();

    // Material goes here
    // material and mesh creates an id for
};



class Scene {
    public:

    // static void LoadMeshes();
    // static void LoadShaders();
    // static void LoadCameras();
    std::unordered_map<std::string, EntityInstanceData> entityInstanceMap;



    std::vector<Entity> entities;

    std::vector<Primitive> debugMeshes;
    std::vector<Shader> shaders;
    std::vector<Camera> initCamConfigs;

    Scene();
    ~Scene();

    Scene(const Scene& other) = delete;
    Scene& operator=(const Scene& other) = delete;

    Scene(Scene&& other);
    Scene& operator=(Scene&& other);


    void DeleteSceneObjects();

    // bool RegisterEntity(Mesh* mesh, Shader* shader);
    // bool RegisterCamera(Camera * camera);

    static Scene GenerateDefaultScene();
    static Scene GenerateBasicScene();

    void RebindAllMeshes();

};

#endif