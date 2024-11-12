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

typedef struct ParticleSystemParameters {
	uint32_t particleCount;
	float timeStep;
	std::string shaderName;

	float fluidViscosity;
	float gravityAccel;
	float maxLifetime;
	float maxSpeed;
};

typedef struct ParticleDataBlock{
		float mass;
		float lifeSpan;
		float radius;
};

typedef struct ParticleSystemDataBlock {
	uint32_t particleCount;
	float timeStep;
	// float fluidViscosity;
	// float gravityAccel;
	// float maxLifetime;
	// float maxSpeed;
	// glm::vec3 initialPosition;
	// glm::vec3 initialVelocity;
};

template <
	typename BaseParticleDataBlock, 
	typename BaseParticleSystemDataBlock, 
	typename BaseParticleSystemParameters
> 
class BaseParticleSystem {
	public: 

	uint32_t particleCount;
	GLuint EBO;
	//std::vector<Vertex> positions;
	GLuint positionsBuffer; // Buffers for freedom
	GLuint velocityBuffer;
	GLuint forcesBuffer;
	// GLuint particleDataBuffer;
	
	std::vector<glm::vec4> positionsVec;
	std::vector<glm::vec4> velocityVec;
	std::vector<glm::vec4> forceVec;

	std::vector<BaseParticleDataBlock> particleDataVec;
	BaseParticleSystemDataBlock particleSystemDataBlock;
	GLuint particleSystemDataBuffer;

	// GLuint particleConstantsDataBuffer;
	Shader* particleShader;
	ComputeShader* particleComputeShader;
	
	BaseParticleSystem(BaseParticleSystemParameters params);
	void InitializeBuffers();
	void CalculateForce() const;

	// void Render();
};


 typedef BaseParticleSystem<
 			ParticleDataBlock, 
 			ParticleSystemDataBlock, 
 			ParticleSystemParameters
 		> ParticleSystem;
#endif