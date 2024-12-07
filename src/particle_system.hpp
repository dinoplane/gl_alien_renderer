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
- TODO Put Perforce and other skills in resume 
- TODO Buy Jo's mp3
*/

 struct ParticleSystemParameters {
	uint32_t particleCount;
	float timeStep;
	std::string shaderName;
};

 struct ParticleDataBlock{
		float mass;
};

 struct ParticleSystemDataBlock {
	uint32_t particleCount;
	float timeStep;
};

class IBaseParticleSystem {
public:
	uint32_t particleCount;
	uint32_t indiceCount;
	GLuint EBO;
	GLuint positionBuffer; // Buffers for freedom
	std::string type;

	virtual void Initialize(void* params) {};
	virtual void InitializeSystemData(void* params) {};
	virtual void InitializeBufferData(void* params) {};
	virtual void InitializeShaders(void* params) {};
	virtual void InitializeBuffers() {};
	virtual void BindBuffers() const{};
	virtual void SetupRender() {};
	virtual void CalculateForces() {};
	virtual void RenderDebug(GLuint VAO) {};
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
	
	std::vector<float> positionVec;
	std::vector<float> velocityVec;
	std::vector<float> forceVec;
    uint32_t vecValCount;

	std::vector<uint32_t> indicesVec;

	std::vector<BaseParticleDataBlock> particleDataVec;
	BaseParticleSystemDataBlock particleSystemDataBlock;
	GLuint particleSystemDataBuffer;

	Shader* particleShader;
	ComputeShader* particleComputeShader;

	virtual void Initialize(void* params) override;
	virtual void InitializeSystemData(void* params) override;
	virtual void InitializeBufferData(void* params) override;
	virtual void InitializeShaders(void* params) override;
	virtual void InitializeBuffers();
	virtual void BindBuffers() const override;
	virtual void CalculateForces() override;
	virtual void RenderDebug(GLuint VAO) override {};
	virtual void SetupRender() override;

};


template class BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters>;
typedef BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters> ParticleSystem;


#endif