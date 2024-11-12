#include <particle_system.hpp>
#include <glm/gtc/random.hpp>

template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> BaseParticleSystem<
	BaseParticleDataBlock, 
	BaseParticleSystemDataBlock, 
	BaseParticleSystemParameters
>::BaseParticleSystem(BaseParticleSystemParameters params){
    // particleCount = params.particleCount;
    // timeStep = params.timeStep;
    // fluidViscosity = params.fluidViscosity;
    // gravityAccel = params.gravityAccel;
    // maxLifetime = params.maxLifetime;
    // maxSpeed = params.maxSpeed;

    particleSystemDataBlock.particleCount = params.particleCount;
    particleSystemDataBlock.timeStep = params.timeStep;
    
    particleCount = params.particleCount;

    positionsVec.resize(params.particleCount);
    velocityVec.resize(params.particleCount);
    forceVec.resize(params.particleCount);
    // particleDataVec.resize(particleCount);

    for (uint32_t i = 0; i < particleCount; ++i){
        positionsVec[i] = glm::vec4(glm::sphericalRand(1.0f), 1.0f);// glm::vec3(0.0f);
        velocityVec[i] = glm::vec4(glm::sphericalRand(1.0f), 1.0f);
        forceVec[i] = glm::vec4(-1.0f);
        // particleDataVec[i].mass = 1.0f;
        // particleDataVec[i].lifeSpan = 1.0f;
        // particleDataVec[i].radius = 1.0f;
    }

    std::string vertPath = "./resources/shader/" + params.shaderName + ".vert";
    std::string fragPath = "./resources/shader/" + params.shaderName + ".frag";
    std::string compPath = "./resources/shader/" + params.shaderName + ".comp";

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
    glCreateBuffers(1, &positionsBuffer);
    glNamedBufferStorage(
        positionsBuffer, 
        sizeof(glm::vec4) * particleCount, 
        positionsVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &velocityBuffer);
    glNamedBufferStorage(
        velocityBuffer, 
        sizeof(glm::vec4) * particleCount, 
        velocityVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );

    glCreateBuffers(1, &forcesBuffer);
    glNamedBufferStorage(
        forcesBuffer, 
        sizeof(glm::vec4) * particleCount, 
        forceVec.data(),
         GL_DYNAMIC_STORAGE_BIT
    );
    
    std::vector<uint32_t> indices(particleCount);
    for (uint32_t i = 0; i < particleCount; ++i){
        indices[i] = i;
    }

    glCreateBuffers(1, &EBO);
    glNamedBufferStorage(
        EBO, 
        sizeof(uint32_t) * particleCount, 
        indices.data(),
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
>::CalculateForce() const{
    
    // Calculate forces
    // for (uint32_t i = 0; i < particleCount; ++i){
    //     forceVec[i] = glm::vec3(0.0f);
    //     for (uint32_t j = 0; j < particleCount; ++j){
    //         if (i != j){
    //             glm::vec3 r = positionsVec[j] - positionsVec[i];
    //             float rLen = glm::length(r);
    //             forceVec[i] += glm::normalize(r) * (particleDataVec[i].mass * particleDataVec[j].mass) / (rLen * rLen);
    //         }
    //     }
    // }
}

template class BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters>;

typedef BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters> ParticleSystem;