#ifndef CLOTH_SYSTEM_H
#define CLOTH_SYSTEM_H

#include <util.h>
#include <particle_system.hpp>
#include <particle_system.cpp>
#include <chrono>
#include <Eigen/SparseCore>

using SparseEntries = std::vector<Eigen::Triplet<double>> ;


#define PARDISO_SOLVE 0

struct ClothSystemParameters {
    uint32_t particleCount;
    double timeStep;
    std::string shaderName;
    uint32_t clothSideLength;
    double cellSideLength;
    double totalMass;
    double gravityAccel;
    double youngModulus;
    double thickness;
};

struct ClothDataBlock{
    double mass;
};

struct ClothSystemDataBlock {
    uint32_t particleCount;
    double timeStep;
    double tolerance;
    double gravityAccel;
    double bendingStiffness;
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
        Eigen::SparseMatrix<double, 
        #if PARDISO_SOLVE == 1
        Eigen::RowMajor
        #else
        Eigen::ColMajor
        #endif
        > massMatrix;
        std::vector<uint32_t> fixedNodes;
        std::vector<uint32_t> fixedIndices;
        std::vector<uint32_t> freeIndices;
        std::vector<uint32_t> oldToNewIndiceMapping;


        struct PickedClothData {
            uint32_t pickedEdgeIdx = 0; 
            uint32_t pickedHingeIdx = 0;
        } pickedClothData;
        GLuint pickedDebugDataBuffer; // picked databuffer
        
        GLuint fixedNodesVAO;
        GLuint fixedNodesDataBuffer;
        uint32_t fixedNodesCount;

        glm::vec4 edgeDebugDataVec[4] { glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f), glm::vec4(0.0f) };
        uint32_t edgeDebugEBOData[6] { 0, 1, 2, 0, 3, 1};
        GLuint edgeDebugDataBuffer; // 4 vec4s
        GLuint edgeDebugEBO;

        Shader* edgeDebugShader;
        Shader* hingeDebugShader;
        Shader* fixedNodesDebugShader;

        bool timerRunning = true;
        double totalForceTime;
        double totalSolveTime;
        
        struct PardisoData {
            void* pt[64];
            int mtype;
            int solver;
            int iparm[65];
            double dparm[64];
            int maxfct, mnum, phase, error, msglvl, num_procs;
            int n , nrhs;
            int perm;


        } pardisoData;

    virtual void InitializeSystemData(void* params) override;
    virtual void InitializeBufferData(void* params) override;
    virtual void InitializeShaders(void* params) override;
    virtual void InitializeBuffers() override;

    virtual void Clear() override;
    virtual void Reinitialize(void* params) override;

    virtual void BindBuffers() const override;
    virtual void SetupRender() override;
    virtual void RenderDebug(GLuint VAO) override;
    virtual void CalculateForces() override;
    void InitializePardiso();
    

    std::chrono::time_point<std::chrono::high_resolution_clock> startTime;
    std::chrono::time_point<std::chrono::high_resolution_clock> lastTime;
    double totalTime;
};


#endif