#include <particle_system.hpp>
#include <glm/gtc/random.hpp>


//template <
//	typename BaseParticleDataBlock, 
//	typename BaseParticleSystemDataBlock, 
//	typename BaseParticleSystemParameters
//> BaseParticleSystem<
//	BaseParticleDataBlock, 
//	BaseParticleSystemDataBlock, 
//	BaseParticleSystemParameters
//>::BaseParticleSystem(BaseParticleSystemParameters* params){
//    // BaseParticleSystemParameters* baseParams = static_cast<BaseParticleSystemParameters*>(params);
//    InitializeSystemData(params); // validate it!
//    InitializeBufferData(params);
//    InitializeShaders(params);
//}

template <
    typename BaseParticleDataBlock, 
    typename BaseParticleSystemDataBlock, 
    typename BaseParticleSystemParameters
> 
void BaseParticleSystem<
    BaseParticleDataBlock, 
    BaseParticleSystemDataBlock, 
    BaseParticleSystemParameters
>::Initialize(void* params){
    BaseParticleSystemParameters* baseParams = static_cast<BaseParticleSystemParameters*>(params);
    InitializeSystemData(baseParams); // validate it!
    InitializeBufferData(baseParams);
    InitializeShaders(baseParams);
}

template <
    typename BaseParticleDataBlock,
    typename BaseParticleSystemDataBlock,
    typename BaseParticleSystemParameters
>
void BaseParticleSystem<
    BaseParticleDataBlock,
    BaseParticleSystemDataBlock,
    BaseParticleSystemParameters
>::InitializeSystemData(void* params) {
    BaseParticleSystemParameters* baseParams = static_cast<BaseParticleSystemParameters*>(params);
    particleSystemDataBlock.particleCount = baseParams->particleCount;
    particleSystemDataBlock.timeStep = baseParams->timeStep;
    particleCount = baseParams->particleCount;
    vecValCount = particleCount * 4;
    indiceCount = particleCount;

}


template <
    typename BaseParticleDataBlock,
    typename BaseParticleSystemDataBlock,
    typename BaseParticleSystemParameters
>
void BaseParticleSystem<
    BaseParticleDataBlock,
    BaseParticleSystemDataBlock,
    BaseParticleSystemParameters
>::InitializeBufferData(void* params) {
    BaseParticleSystemParameters* baseParams = static_cast<BaseParticleSystemParameters*>(params);
    positionVec.resize(vecValCount);
    velocityVec.resize(vecValCount);
    forceVec.resize(vecValCount);
    // particleDataVec.resize(particleCount);
    fmt::printf("KILLAKILL");

    for (uint32_t i = 0; i < particleCount; ++i) {
        glm::vec3 tmpPosition = glm::sphericalRand(1.0f);
        glm::vec3 tmpVelocity = glm::sphericalRand(1.0f);

        for (uint32_t j = 0; j < 3; ++j) {
            positionVec[i * 4 + j] = tmpPosition[j];
            velocityVec[i * 4 + j] = tmpVelocity[j];
            forceVec[i * 4 + j] = -1.0f;
        }
        positionVec[i * 4 + 3] = 1.0f;
        velocityVec[i * 4 + 3] = 1.0f;
        forceVec[i * 4 + 3] = 1.0f;

        // particleDataVec[i].mass = 1.0f;
        // particleDataVec[i].lifeSpan = 1.0f;
        // particleDataVec[i].radius = 1.0f;
    }

    fmt::print("InistializeBufferData");
    indicesVec.resize(indiceCount);
    for (uint32_t i = 0; i < indiceCount; ++i) {
        indicesVec[i] = i;
    }
    

}


template <
    typename BaseParticleDataBlock,
    typename BaseParticleSystemDataBlock,
    typename BaseParticleSystemParameters
>
void BaseParticleSystem<
    BaseParticleDataBlock,
    BaseParticleSystemDataBlock,
    BaseParticleSystemParameters
>::InitializeShaders(void* params) {
    BaseParticleSystemParameters* baseParams = static_cast<BaseParticleSystemParameters*>(params);
    std::string vertPath = "./resources/shader/" + baseParams->shaderName + "_particle.vert";
    std::string fragPath = "./resources/shader/" + baseParams->shaderName + "_particle.frag";
    std::string compPath = "./resources/shader/" + baseParams->shaderName + "_particle.comp";

    particleShader = new Shader(vertPath.c_str(), fragPath.c_str());
    particleComputeShader = new ComputeShader(compPath.c_str());
}


template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> 
void BaseParticleSystem<
	BaseParticleDataBlock, 
	BaseParticleSystemDataBlock, 
	BaseParticleSystemParameters
>::InitializeBuffers(){
    glCreateBuffers(1, &positionBuffer);
    glNamedBufferStorage(
        positionBuffer, 
        sizeof(float) * vecValCount, 
        positionVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &velocityBuffer);
    glNamedBufferStorage(
        velocityBuffer, 
        sizeof(float) * vecValCount, 
        velocityVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &forcesBuffer);
    glNamedBufferStorage(
        forcesBuffer, 
        sizeof(float) * vecValCount, 
        forceVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );
    

    glCreateBuffers(1, &EBO);
    glNamedBufferStorage(
        EBO, 
        sizeof(uint32_t) * indiceCount, 
        indicesVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    // glCreateBuffers(1, &particleDataBuffer);
    // glNamedBufferStorage(
    //     particleDataBuffer, 
    //     sizeof(ParticleDataBlock) * particleCount, 
    //     nullptr,
    //      GL_DYNAMIC_STORAGE_BIT
    // );

     glCreateBuffers(1, &particleSystemDataBuffer);
     glNamedBufferStorage(
         particleSystemDataBuffer, 
         sizeof(BaseParticleSystemDataBlock), 
         &particleSystemDataBlock,
          GL_DYNAMIC_STORAGE_BIT
     );
}

template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> 
void BaseParticleSystem<
	BaseParticleDataBlock, 
	BaseParticleSystemDataBlock, 
	BaseParticleSystemParameters
>::CalculateForces() {
    
    //BindParticleSystem(&particleSystem);
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

template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> 
void BaseParticleSystem<
	BaseParticleDataBlock, 
	BaseParticleSystemDataBlock, 
	BaseParticleSystemParameters
>::BindBuffers() const {
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_POSITIONS_SSBO_BINDING, positionBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_VELOCITIES_SSBO_BINDING, velocityBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FORCES_SSBO_BINDING, forcesBuffer);
    //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_DATA_SSBO_BINDING, particleDataBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_SYSTEM_DATA_SSBO_BINDING, particleSystemDataBuffer);
    particleComputeShader->use();
}

template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> 
void BaseParticleSystem<
	BaseParticleDataBlock, 
	BaseParticleSystemDataBlock, 
	BaseParticleSystemParameters
>::SetupRender(){

    particleShader->use();

}


