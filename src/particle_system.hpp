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

DONE? redefine the particle system design
DONE! Get the cloth mesh to render / show up on the screen
DONE! set up hinges correctly
DONE! set up edges correctly
DONE! calculate the gradient of stretching
DONE! calculate the hessian of stretching
DONE! calculate the gradient of bending (hinges)
	- DONE! calculate gradTheta
	- DONE! calculate getTheta
	- DONE! Debug the gradTheta
DONE! calculate the hessian of bending (hinges)
	- DONE! calculate hessTheta

- DONE! Debug the gradient 
- DONE! Debug the hessian

- DONE! get the cloth mesh to update

- DONE! translate the mainloop
	DONE! translate the stretching calculation
	- DONE! translate the bending calculation
	- DONE! translate the rest of the calculation (the solver)
	- DONE! mess with free indices
	- DONE! translate the error loop

- CUT create compute shaders for the operations
- CUT use a gpu to solve these

- DONE Record Current Times
- DONE Profile which areas need optimization
    - TODO Find the hot paths
	- TODO find the solutions
	- TODO implement the solutions

- TODO Change the final hessians into Sparse Matrices
	- TODO Stretching
	- TODO Bending
	- TODO Final
	- TODO Write an algorithm to offset indices with respect to the free indices

- TODO Integrate the pardiso solver
- TODO Change the Edge and Hinge layouts to prefer a column major layout
- TODO Multithread the gradient and hessian calculations
	- Fork threads
	- Thread pool!
- TODO Use Matrices when when possible...


- TODO create a presentation
	- DONE? midterm (3 slides?)
	- TODO final (5 slides)
- TODO create a report
	- DONE midterm (5 pages, 5+ references)
	- TODO final (n pages, 10 references)


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
	GLuint positionBuffer; // Buffers for freedom
	std::string type;

	// IBaseParticleSystem(void* params) {};
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
	// GLuint particleDataBuffer;
	
	std::vector<float> positionVec;
	std::vector<float> velocityVec;
	std::vector<float> forceVec;
    uint32_t vecValCount;

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
	virtual void CalculateForces() override;
	virtual void RenderDebug(GLuint VAO) override {};
	virtual void SetupRender() override;

	// void Render();
};


template class BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters>;
typedef BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters> ParticleSystem;


#endif