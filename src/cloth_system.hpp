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
        uint32_t freeDOFCount;

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
        GLuint fixedNodesVAO;
        GLuint fixedNodesDataBuffer;
        uint32_t fixedNodesCount;


        // GLuint fixedNodesEBO;

        glm::vec4 edgeDebugDataVec[4] { glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f) };
        uint32_t edgeDebugEBOData[6] { 0, 1, 2, 0, 3, 1};
        GLuint edgeDebugDataBuffer; // 4 vec4s
        GLuint edgeDebugEBO;

        Shader* edgeDebugShader;
        Shader* hingeDebugShader;
        Shader* fixedNodesDebugShader;

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
    
/*
Timing Data:

Debug Build

c    Single threaded             Edges  Hinges
1    4 nodes  : 0.1 ms           5      1
2    9 nodes  : 0.9 ms          16      8
3    16 nodes : 24-47 ms        33      21
4    25 nodes : 45-61 ms        56      40
5    36 Nodes : 70-134 ms       85      65
6    49 nodes : 160-227 ms      120     96

Release build (O2 optimization)
c    Single threaded      
1    4 nodes   : 0.1 ms    
2    9 nodes   : 0.2 ms    
3    16 nodes  : 0.9 ms  
4    25 nodes  : 1.04 ms  
5    36 Nodes  : 1.05 ms 
6    49 nodes  : 0.84 ms
7    64 nodes  : 0.827 ms
8    81 nodes  : 1.338 ms
9    100 nodes : 1.449 ms
10   121 nodes : 11.96 ms
11   144 nodes : 33.754 ms
12   169 nodes : 46.387 ms
13   196 nodes : 64.464 ms
14   225 nodes : 89.895  ms
15   256 nodes : 118.725 ms
16   289 nodes : 162.788 ms
17   324 nodes : 214.5 ms
18   361 nodes : 346.022 ms




*/


#endif