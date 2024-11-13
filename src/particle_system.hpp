#ifndef PARTICLE_SYSTEM_H
#define PARTICLE_SYSTEM_H

#include <util.h>

#include<fastgltf/core.hpp>
#include <mesh.hpp>
#include <transform.hpp>
#include <volume.hpp>
#include <vector>
#include <glm/glm.hpp>
#include <shader_s.hpp>
#include <shader_c.hpp>
#include <gl_bindings.h>

/* TODOs
- redefine the particle system design
- Get the cloth mesh to render / show up on the screen
- get the cloth mesh to update
- set up hinges correctly
- set up edges correctly
- calculate the gradient of theta (hinges)
- calculate the hessian of theta (hinges)
- calculate the gradient of stretching
- calculate the hessian of stretching
- create compute shaders for the operations
- use a gpu to solve these
- create a presentation
- create a report


*/

 struct ParticleSystemParameters {
	uint32_t particleCount;
	float timeStep;
	std::string shaderName;

	// float fluidViscosity;
	// float gravityAccel;
	// float maxLifetime;
	// float maxSpeed;
};

 struct ParticleDataBlock{
		float mass;
		// float lifeSpan;
		// float radius;
};

 struct ParticleSystemDataBlock {
	uint32_t particleCount;
	float timeStep;
	// float fluidViscosity;
	// float gravityAccel;
	// float maxLifetime;
	// float maxSpeed;
	// glm::vec3 initialPosition;
	// glm::vec3 initialVelocity;
};

class IBaseParticleSystem {
public:
	uint32_t particleCount;
	uint32_t indiceCount;
	GLuint EBO;
	//std::vector<Vertex> positions;
	GLuint positionsBuffer; // Buffers for freedom

	// IBaseParticleSystem(void* params) {};
	virtual void Initialize(void* params) {};
	virtual void InitializeSystemData(void* params) {};
	virtual void InitializeBufferData(void* params) {};
	virtual void InitializeShaders(void* params) {};
	virtual void InitializeBuffers() {};
	virtual void BindBuffers() const{};
	virtual void SetupRender() {};
	virtual void CalculateForces() const {};
	virtual ~IBaseParticleSystem() = default;
};

template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> 
class BaseParticleSystem : public IBaseParticleSystem{
	public: 


	GLuint velocityBuffer;
	GLuint forcesBuffer;
	// GLuint particleDataBuffer;
	
	std::vector<glm::vec4> positionsVec;
	std::vector<glm::vec4> velocityVec;
	std::vector<glm::vec4> forceVec;
	std::vector<uint32_t> indicesVec;

	std::vector<BaseParticleDataBlock> particleDataVec;
	BaseParticleSystemDataBlock particleSystemDataBlock;
	GLuint particleSystemDataBuffer;

	// GLuint particleConstantsDataBuffer;
	Shader* particleShader;
	ComputeShader* particleComputeShader;
	//BaseParticleSystem();
	//BaseParticleSystem();
	virtual void Initialize(void* params) override;
	virtual void InitializeSystemData(void* params) override;
	virtual void InitializeBufferData(void* params) override;
	virtual void InitializeShaders(void* params) override;
	virtual void InitializeBuffers();
	virtual void BindBuffers() const override;
	virtual void CalculateForces() const override;
	virtual void SetupRender() override;

	// void Render();
};


template class BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters>;
typedef BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters> ParticleSystem;


#endif