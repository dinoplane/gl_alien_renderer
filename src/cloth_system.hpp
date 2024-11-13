#ifndef CLOTH_SYSTEM_H
#define CLOTH_SYSTEM_H

#include <util.h>
// #include <fastgltf/core.hpp>
#include <particle_system.hpp>
#include <particle_system.cpp>

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

}; // Add as we figure that out

struct EdgeData{
 
}; // Add as we figure this out



class ClothSystem : public BaseParticleSystem<
    ClothDataBlock,
    ClothSystemDataBlock,
    ClothSystemParameters
> {
    public:

        std::vector<glm::vec4> normalsVec;
        GLuint normalsBuffer;
        std::vector<glm::vec4> edgesVec;
        GLuint edgesBuffer;
        std::vector<glm::vec4> hingesVec;
        GLuint hingesBuffer;

    //ClothSystem(ClothSystemParameters* params);
        //virtual void Initialize(void* params) override;
    virtual void InitializeSystemData(void* params) override;
    virtual void InitializeBufferData(void* params) override;
    virtual void InitializeShaders(void* params) override;
    virtual void InitializeBuffers() override;
    virtual void BindBuffers() const override;
    virtual void SetupRender() override;
    virtual void CalculateForces() const override;
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