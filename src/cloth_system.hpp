#ifndef CLOTH_SYSTEM_H
#define CLOTH_SYSTEM_H

#include <util.h>
// #include <fastgltf/core.hpp>
#include <particle_system.hpp>
#include <particle_system.cpp>
#include <chrono>

struct ClothSystemParameters {
    uint32_t particleCount;
    float timeStep;
    std::string shaderName;
     uint32_t clothSideLength;
     float cellSideLength;
    //float tolerance;
    float totalMass;
    //float fluidViscosity;
    //float gravityAccel;
    //float initialPosition;

};

struct ClothDataBlock{
    float mass;
    // bool fixedIndex; // dunno if that should be here...
};

struct ClothSystemDataBlock {
    uint32_t particleCount;
    float timeStep;
    // float tolerance;
    // float fluidViscosity;
    // float gravityAccel;
    // float youngModulus;
    // float thickness;
    // float bendingStiffness;
    // float stretchStiffness;
    
};

template class BaseParticleSystem<ClothDataBlock, ClothSystemDataBlock, ClothSystemParameters>;



struct HingeData{
    glm::ivec4 nodes;
}; // Add as we figure that out

struct EdgeData{
    glm::ivec2 nodes;
}; // Add as we figure this out



class ClothSystem : public BaseParticleSystem<
    ClothDataBlock,
    ClothSystemDataBlock,
    ClothSystemParameters
> {
    uint32_t dofCount;
    public:

        std::vector<glm::vec4> normalsVec;
        GLuint normalsBuffer;
        
        std::vector<glm::ivec2> edgeVec;
        uint32_t edgeCount;
        GLuint edgesBuffer;

        std::vector<glm::ivec4> hingeVec;
        uint32_t hingeCount;
        GLuint hingesBuffer;


        struct PickedClothData {
            uint32_t pickedEdgeIdx = 0; 
            uint32_t pickedHingeIdx = 0;
        } pickedClothData;
        GLuint pickedDebugDataBuffer; // picked databuffer
        
        // glm::vec4 edgePointDataVec[2] { glm::vec4(0.0f), glm::vec4(0.0f) } ; // picked data vec
        // uint32_t edgeDebugEBOData[2] { 0, 1 };
        // GLuint edgePointDebugBuffer; // 2 vec4s
        // GLuint edgePointEBO;

        glm::vec4 edgeDebugDataVec[4] { glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f) };
        uint32_t edgeDebugEBOData[6] { 0, 1, 2, 0, 3, 1};
        GLuint edgeDebugDataBuffer; // 4 vec4s
        GLuint edgeDebugEBO;

        Shader* edgeDebugShader;
        Shader* hingeDebugShader;

        

    //ClothSystem(ClothSystemParameters* params);
        //virtual void Initialize(void* params) override;
    virtual void InitializeSystemData(void* params) override;
    virtual void InitializeBufferData(void* params) override;
    virtual void InitializeShaders(void* params) override;
    virtual void InitializeBuffers() override;
    virtual void BindBuffers() const override;
    virtual void SetupRender() override;
    virtual void RenderDebug(GLuint VAO) override;
    virtual void CalculateForces() const override;

    std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;

    // void CalculateNewPositions();
};
//class ClothSystem : public BaseParticleSystem<
//    ClothDataBlock, 
//    ClothSystemDataBlock, 
//    ClothSystemParameters
//> {
//    public: 
//

//
//    ClothSystem(ClothSystemParameters* params);
//    void InitializeSystemData(void* params) override;
//    void InitializeBufferData(void* params) override;
//    //void InitializeShaders(void* params) override;
//    //void InitializeBuffers() override;
//    //void BindBuffers() const override;
//    //void SetupRender() override;
//    //void CalculateForces() const override;
//};
    
    
#endif