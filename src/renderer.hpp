#ifndef RENDERER_H
#define RENDERER_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <vector>
#include <gpu_structs.hpp>
// #include <memory>
#include <shader_s.hpp>
#include <shader_c.hpp>

#include <unordered_map>
#include <particle_system.hpp>

class Scene;
class Camera;
struct Primitive;
struct EntityInstanceData;
struct ComputeShader;

// class Light;
// class Material;
// class Entity;

// Struct for MultiDrawElements
typedef struct {
    unsigned int  count;
    unsigned int  instanceCount;
    unsigned int  firstIndex;
    int           baseVertex;
    unsigned int  baseInstance;
    // Optional user-defined data goes here - if nothing, stride is 0
} DrawElementsIndirectCommand;
// sizeof(DrawElementsIndirectCommand) == 20

class Renderer {
    public:

    static std::vector<Camera> allCameras;
    static std::vector<Primitive> allDebugMeshes;
    static Shader* debugShader;
    static Shader* debugWireShader;

    static ComputeShader* cullShader;
    static ComputeShader* instCountSetShader;

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
    static void ToggleCull(const Scene& scene);



    GLuint VAO;
    GLuint debugVAO;
    GLuint instVAO;
    GLuint particleVAO;

    static GLuint quadVAO;
    static GLuint quadVBO;


    // we can make a more elaborate system later when we want to make the canny edge detection
    GLuint FBO;
    GLuint srcFBOTexture; // texture to use dsa
    GLuint RBO;

    // GLuint FBO;
    GLuint dstFBOTexture; // texture to use dsa
    // GLuint RBO;


    GLuint cameraMatricesUBO;
    struct CameraMatricesUBOBlock {
        glm::mat4 projection;
        glm::mat4 view;
    } cameraMatricesUBOBlock;

    GLuint frustumCullDataUBO;
    struct alignas(16) FrustumCullDataUBOBlock {
        GPUFrustum frustum;
        GPUSphere boundingVolume;
        uint32_t instCount;
        uint32_t doCull;
    } frustumCullDataUBOBlock;
    static uint32_t doCull;



    // GLuint meshPropertiesUBO;
    // struct MeshPropertiesUBOBlock {
    //     glm::mat4 modelFromMesh;
    // } meshPropertiesUBOBlock;

    GLuint materialUniformsUBO;
    // struct MaterialUniformsUBOBlock {
    //     glm::vec4 baseColorFactor;
    //     float alphaCutoff;
    //     uint32_t flags;
    // } materialUniformsUBOBlock;

    // 1 draw command per instance
    std::vector<DrawElementsIndirectCommand> commands;

    GLuint inDrawCmdBuffer;
    GLuint outDrawCmdBuffer;


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
    void CreateParticleVAO();
    void CreateDebugVAO();
    void CreateFBO(float w, float h);
    void CreateRBO(float w, float h);

    void Resize(float w, float h);
    // void Update(const Scene& scene);

    // void BindMesh(const Primitive& mesh);
    void BindPrimitive(const Primitive& mesh);

    void BindInstanceData(const EntityInstanceData& entityInstanceData);
    void BindInstancePrimitive(const Primitive& primitive);
    void BindParticleSystem(const IBaseParticleSystem* particleSystem);

    void BindDebugMesh(const Primitive& mesh);

    // void Render(Scene* scene);


    void Render(const Scene& scene);
    void RenderEntities(const Scene& scene);
    void RenderInstancedStaticModels(const Scene& scene);
    void RenderParticleSystems(const Scene& scene);
    void RenderDebugVolumes(const Scene& scene);
    
    void RenderPostProcess();


    ComputeShader* frustumCullShader;
    // frustum is uniform

    // model transform buffer inside scene

    void BindInstanceCullingBuffers(const EntityInstanceData& entityInstanceData);
    void BindDrawCmdInstCountSetBuffers(const EntityInstanceData& entityInstanceData);

    void SwitchCamera(const Scene& scene);

    public:
    uint32_t mainCameraIdx {0};
    uint32_t currCameraIdx {0};
    bool debugOn {true};
};

#endif