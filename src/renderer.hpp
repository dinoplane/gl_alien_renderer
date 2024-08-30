#ifndef RENDERER_H
#define RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
// #include <memory>
#include <shader_s.hpp>
#include <unordered_map>
class Scene;
class Camera;
struct Mesh;
struct EntityInstanceData;
struct ComputeShader;
// class Light;
// class Material;
// class Entity;

class Renderer {
    public:

    static std::vector<Camera> allCameras;
    static std::vector<Mesh> allDebugMeshes;
    static Shader* debugShader;
    static Shader* debugWireShader;



    static std::unordered_map<std::string, Shader*> allPostProcessShaders;
    static Shader* postProcessShader;
    static Shader* passthroughShader;

    static constexpr float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
        // positions   // texCoords
        -1.0f,  1.0f,  0.0f, 1.0f,
        -1.0f, -1.0f,  0.0f, 0.0f,
        1.0f, -1.0f,  1.0f, 0.0f,

        -1.0f,  1.0f,  0.0f, 1.0f,
        1.0f, -1.0f,  1.0f, 0.0f,
        1.0f,  1.0f,  1.0f, 1.0f
    };



    static void SetupScreenQuad();



    GLuint VAO;
    GLuint debugVAO;
    GLuint instVAO;

    static GLuint quadVAO;
    static GLuint quadVBO;


    // we can make a more elaborate system later when we want to make the canny edge detection
    GLuint FBO;
    GLuint srcFBOTexture; // texture to use dsa
    GLuint RBO;

    // GLuint FBO;
    GLuint dstFBOTexture; // texture to use dsa
    // GLuint RBO;


    float width;
    float height;

    // static void LoadMeshes();
    // static void LoadShaders();
    // static void LoadCameras();

    // Scene* scene;
    Renderer(float w, float h);

    void Init(float w, float h);
    void CreateVAO();
    void CreateInstanceVAO();
    void CreateDebugVAO();
    void CreateFBO(float w, float h);
    void CreateRBO(float w, float h);

    void Resize(float w, float h);
    // void Update(const Scene& scene);

    void BindMesh(Mesh* mesh);
    void BindInstanceMesh(EntityInstanceData* entityInstanceData);
    void BindDebugMesh(Mesh* mesh);

    void Render(Scene* scene);
    void RenderPostProcess();


    ComputeShader* frustumCullShader;
    // frustum is uniform

    // model transform buffer inside scene


    void SetupCullStructures();


    void SwitchCamera(const Scene& scene);

    public:
    uint mainCameraIdx {0};
    uint currCameraIdx {0};
    bool debugOn {true};
};

#endif