#ifndef CLOTH_SYSTEM_H
#define CLOTH_SYSTEM_H

#include <util.h>
// #include <fastgltf/core.hpp>
#include <particle_system.hpp>
#include <particle_system.cpp>
#include <chrono>

struct ClothSystemParameters {
    uint32_t particleCount;
    double timeStep;
    std::string shaderName;
     uint32_t clothSideLength;
    double cellSideLength;
    double totalMass;
//  double fluidViscosity;
    double gravityAccel;
//  double initialPosition;
    double youngModulus;
    double thickness;


};

struct ClothDataBlock{
    double mass;
    // bool fixedIndex; // dunno if that should be here...
};

struct ClothSystemDataBlock {
    uint32_t particleCount;
    double timeStep;
    double tolerance;
    double gravityAccel;
    double bendingStiffness;
    // double fluidViscosity;
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
    public:
        std::vector<float> clothPositionVec;
        std::vector<double> clothVelocityVec;
        std::vector<double> clothForceVec;


        uint32_t dofCount;

        std::vector<glm::vec4> normalsVec;
        GLuint normalsBuffer;
        
        std::vector<glm::ivec2> edgeVec;
        std::vector<double> undeformedEdgeLengthVec;
        std::vector<double> elasticStretchingVec;
        double bendingStiffness;
        double dt;
        double gravityAccel;
        double tol;

        uint32_t edgeCount;
        GLuint edgesBuffer;

        std::vector<glm::ivec4> hingeVec;
        uint32_t hingeCount;
        GLuint hingesBuffer;
        std::vector<double> thetaBarVec;

        std::vector<double> externalForcesVec;
        
        std::vector<double> dofPositions;
        std::vector<double> dofVelocities;
        std::vector<double> lastdofPositions;

        std::vector<double> massVector;
        std::vector<double> massMatrix;
        std::vector<uint32_t> fixedNodes;
        std::vector<uint32_t> fixedIndices;
        std::vector<uint32_t> freeIndices;


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
        bool timerRunning = true;
        

    //ClothSystem(ClothSystemParameters* params);
        //virtual void Initialize(void* params) override;
    virtual void InitializeSystemData(void* params) override;
    virtual void InitializeBufferData(void* params) override;
    virtual void InitializeShaders(void* params) override;
    virtual void InitializeBuffers() override;
    virtual void BindBuffers() const override;
    virtual void SetupRender() override;
    virtual void RenderDebug(GLuint VAO) override;
    virtual void CalculateForces() override;

    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
    double totalTime;


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