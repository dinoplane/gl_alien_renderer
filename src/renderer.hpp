#ifndef RENDERER_H
#define RENDERER_H

#include <vector>
// #include <memory>

class Scene;
class Camera;
// class Light;
// class Material;
// class Entity;

class Renderer {
    public:

    // static void LoadMeshes();
    // static void LoadShaders();
    // static void LoadCameras();

    // Scene* scene;
    void Init();
    void Resize(const Camera& camera);
    // void Update(const Scene& scene);
    void Render(Scene* scene);

    void SwitchCamera(const Scene& scene);

    public:
    uint mainCameraIdx {0};
    bool debugOn {true};
};

#endif