#ifndef SCENE_H
#define SCENE_H

#include <vector>
#include <unordered_map>
// #include <memory>
#include <mesh.hpp>
class Entity;
class Shader; // Change this to only use move semantics
class Camera;
struct GPUSphere;
// class Light;
// class Material;
// class Entity;

struct EntityInstanceData
{
    Mesh instMesh;
    uint instCount;
    std::vector<glm::mat4> modelToWorldMat;
    GLuint instModelMatrixBuffer;
    std::vector<GPUSphere> boundingVolumes;
    GLuint instBoundingVolumeBuffer;
    GLuint instMeshRenderedBuffer;


    void GenerateInstanceBuffers();

    // Material goes here
    // material and mesh creates an id for
};

class Scene {
    public:

    // static void LoadMeshes();
    // static void LoadShaders();
    // static void LoadCameras();
    std::unordered_map<uint, EntityInstanceData> entityInstanceMap;



    std::vector<Entity> entities;

    std::vector<Mesh> debugMeshes;
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