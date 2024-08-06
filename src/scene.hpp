#ifndef SCENE_H
#define SCENE_H

#include <vector>
// #include <memory>
struct Mesh;
class Entity;
class Shader; // Change this to only use move semantics
class Camera;
// class Light;
// class Material;
// class Entity;

class Scene {
    public:

    // static void LoadMeshes();
    // static void LoadShaders();
    // static void LoadCameras();

    std::vector<Entity> entities;
    std::vector<Mesh> debugMeshes;
    std::vector<Shader> shaders;
    std::vector<Camera> cameras;

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

    void RebindAllMeshes();

};

#endif