#include <cloth_system.hpp>
#include <shader_c.hpp>
#include <shader_s.hpp>


//ClothSystem::ClothSystem(ClothSystemParameters* params) :
//    BaseParticleSystem(params) {
//
//};

void ClothSystem::InitializeSystemData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    clothParams->particleCount = (clothParams->clothSideLength + 1) * (clothParams->clothSideLength + 1);
    particleSystemDataBlock.particleCount = clothParams->particleCount;
    particleCount = clothParams->particleCount;
    particleSystemDataBlock.timeStep = clothParams->timeStep;
    //particleSystemDataBlock.tolerance = clothParams->tolerance;
    //particleSystemDataBlock.fluidViscosity = clothParams->fluidViscosity;
    //particleSystemDataBlock.gravityAccel = clothParams->gravityAccel;
    //particleSystemDataBlock.youngModulus = 1.0e7f;
    //particleSystemDataBlock.thickness = 1.0e3f;
    //particleSystemDataBlock.bendingStiffness = 2.0f / sqrt(3.0f) * particleSystemDataBlock.youngModulus * pow(particleSystemDataBlock.thickness, 3.0f) / 12.0f;
    //particleSystemDataBlock.stretchStiffness = 1.0f;
};



void ClothSystem::InitializeBufferData(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
        clothParams->particleCount = (clothParams->clothSideLength+ 1) * (clothParams->clothSideLength + 1);
    // count particles
    positionsVec.resize(clothParams->particleCount);
    velocityVec.resize(clothParams->particleCount);
    forceVec.resize(clothParams->particleCount);
    normalsVec.resize(clothParams->particleCount);
    particleDataVec.resize(clothParams->particleCount);

    // tangentsVec.resize(clothParams->particleCount);
    // bitangentsVec.resize(clothParams->particleCount);
    edgesVec.resize(clothParams->particleCount);
    hingesVec.resize(clothParams->particleCount);

    glm::vec4 meshStartPosition = glm::vec4(0.0f, 0.0f, 0.0f, 1.0f);
    float yIncrement = 0.5f * sqrt(3.0f) * clothParams->cellSideLength;
    uint32_t clothSideParticleCount = clothParams->clothSideLength + 1;
    float deltaMass = clothParams->totalMass / clothParams->particleCount;

    for (uint32_t rowIdx = 0; rowIdx < clothSideParticleCount; ++rowIdx) {
        glm::vec4 rowStartPosition =
            meshStartPosition;

        rowStartPosition.z -= yIncrement * rowIdx;
        if (rowIdx % 2 == 1) {
            rowStartPosition.x += 0.5f * clothParams->cellSideLength;
        }

        glm::vec4 currentPosition = rowStartPosition;
        for (uint32_t colIdx = 0; colIdx < clothSideParticleCount; ++colIdx) {
            uint32_t particleIdx = rowIdx * clothSideParticleCount + colIdx;
            positionsVec[particleIdx] = currentPosition;
            velocityVec[particleIdx] = glm::vec4(0.0f);
            forceVec[particleIdx] = glm::vec4(0.0f);
            particleDataVec[particleIdx].mass = deltaMass;

            // normalsVec[particleIdx] = glm::vec4(0.0f);
            currentPosition.x += clothParams->cellSideLength;
        }
    }

    // Fill in the EBO
    indicesVec.resize(6 * clothParams->clothSideLength * clothParams->clothSideLength);
    for (uint32_t rowIdx = 0; rowIdx < clothParams->clothSideLength; ++rowIdx) {
        for (uint32_t colIdx = 0; colIdx < clothParams->clothSideLength; ++colIdx) {
            // ACD / ABD
            //    C---------D
            //     \   ____/ \
            //      \ /       \
            //       A---------B
            //
            //    ABC / BDC
            //       C---------D
            //      / \____   /
            //     /       \ /
            //    A---------B

            

            uint32_t cIdx = rowIdx * clothSideParticleCount + colIdx;
            uint32_t dIdx = cIdx + 1;
            uint32_t aIdx = (rowIdx + 1) * clothSideParticleCount + colIdx;
            uint32_t bIdx = aIdx + 1;

            uint32_t indiceIdx = 6 * (rowIdx * clothParams->clothSideLength + colIdx);
            fmt::print("\ni {}\n", indiceIdx);
            if (rowIdx % 2 == 0) {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = dIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = aIdx;
                indicesVec[indiceIdx + 4] = bIdx;
                indicesVec[indiceIdx + 5] = dIdx;
            }
            else {
                indicesVec[indiceIdx + 0] = aIdx;
                indicesVec[indiceIdx + 1] = bIdx;
                indicesVec[indiceIdx + 2] = cIdx;

                indicesVec[indiceIdx + 3] = bIdx;
                indicesVec[indiceIdx + 4] = dIdx;
                indicesVec[indiceIdx + 5] = cIdx;
            }

            fmt::print("\nFirst Triangle\n");
            for (uint32_t j = 0; j < 3; ++j) {
                fmt::print("{} ", indicesVec[indiceIdx + j]);
            }
            fmt::print("\nSecond Triangle\n");
            for (uint32_t j = 3; j < 6; ++j) {
                fmt::print("{} ", indicesVec[indiceIdx + j]);
            }
        }
    }
    indiceCount = 6 * clothParams->clothSideLength * clothParams->clothSideLength;
}

void ClothSystem::InitializeShaders(void* params) {
    ClothSystemParameters* clothParams = static_cast<ClothSystemParameters*>(params);
    std::string vertPath = "./resources/shader/" + clothParams->shaderName + ".vert";
    std::string fragPath = "./resources/shader/" + clothParams->shaderName + ".frag";
    std::string compPath = "./resources/shader/" + clothParams->shaderName + ".comp";

    particleShader = new Shader(vertPath.c_str(), fragPath.c_str());
    particleComputeShader = new ComputeShader(compPath.c_str());
}

void ClothSystem::InitializeBuffers(){
    BaseParticleSystem::InitializeBuffers();
    //glCreateBuffers(1, &edgesBuffer);
    //glNamedBufferStorage(
    //    edgesBuffer, 
    //    sizeof(glm::vec4) * particleCount, 
    //    edgesVec.data(),
    //     GL_DYNAMIC_STORAGE_BIT
    //);

    //glCreateBuffers(1, &hingesBuffer);
    //glNamedBufferStorage(
    //    hingesBuffer, 
    //    sizeof(glm::vec4) * particleCount, 
    //    hingesVec.data(),
    //     GL_DYNAMIC_STORAGE_BIT
    //);
}

void ClothSystem::CalculateForces() const{
    // Iterate through every hinge

    // Iterate through every edge

    glNamedBufferSubData(
        particleSystemDataBuffer,
        0,
        sizeof(ParticleSystemDataBlock),
        &particleSystemDataBlock
    );
    glDispatchCompute(particleCount, 1u, 1u);
    glMemoryBarrier(GL_ALL_BARRIER_BITS);

    glFinish();


}

void ClothSystem::BindBuffers() const{
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_POSITIONS_SSBO_BINDING, positionsBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_VELOCITIES_SSBO_BINDING, velocityBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FORCES_SSBO_BINDING, forcesBuffer);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_DATA_SSBO_BINDING, particleDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_SYSTEM_DATA_SSBO_BINDING, particleSystemDataBuffer);
    particleComputeShader->use();
}

void ClothSystem::SetupRender() {

    particleShader->use();
}

