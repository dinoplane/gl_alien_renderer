#ifndef RENDERER_H
#define RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
// #include <memory>
#include <shader_s.hpp>

class Scene;
class Camera;
struct Mesh;
// class Light;
// class Material;
// class Entity;

class Renderer {
    public:

    static std::vector<Camera> allCameras;
    static std::vector<Mesh> allDebugMeshes;
    static Shader* debugShader;
    static Shader* debugWireShader;


    GLuint VAO;
    GLuint debugVAO;

    GLuint FBO;
    GLuint FBOTexture; // texture to use dsa
    GLuint RBO;

    float width;
    float height;

    // static void LoadMeshes();
    // static void LoadShaders();
    // static void LoadCameras();

    // Scene* scene;
    Renderer(float w, float h);

    void Init(float w, float h);
    void CreateVAO();
    void CreateDebugVAO();
    void CreateFBO(float w, float h);
    void CreateRBO(float w, float h);

    void Resize(float w, float h);
    // void Update(const Scene& scene);

    void BindMesh(Mesh* mesh);
    void BindDebugMesh(Mesh* mesh);

    void Render(Scene* scene);

    void SwitchCamera(const Scene& scene);

    public:
    uint mainCameraIdx {0};
    uint currCameraIdx {0};
    bool debugOn {true};
};

#endif