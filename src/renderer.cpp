#include <renderer.hpp>
#include <camera.hpp>
#include <scene.hpp>
#include <shader_s.hpp>
#include <mesh.hpp>
#include <material.hpp>
#include <shader_c.hpp>
#include <entity.hpp>
#include <util.h>
#include <particle_system.hpp>
#include <gl_bindings.h>


#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>

std::vector<Camera> Renderer::allCameras;
std::vector<Primitive> Renderer::allDebugMeshes;
Shader* Renderer::debugShader = nullptr;
Shader* Renderer::debugWireShader = nullptr;
Shader* Renderer::postProcessShader = nullptr;
Shader* Renderer::passthroughShader = nullptr;
ComputeShader* Renderer::cullShader = nullptr;
ComputeShader* Renderer::instCountSetShader = nullptr;

GLuint Renderer::quadVAO = 0;
GLuint Renderer::quadVBO = 0;
uint32_t Renderer::doCull = 0;


void Renderer::SetupScreenQuad(){
    glCreateVertexArrays(1, &Renderer::quadVAO);
    glEnableVertexArrayAttrib(Renderer::quadVAO, 0);
    glEnableVertexArrayAttrib(Renderer::quadVAO, 1);

    // position attribute
    glVertexArrayAttribFormat(Renderer::quadVAO, 0, 2, GL_FLOAT, GL_FALSE, 0);
    glVertexArrayAttribBinding(Renderer::quadVAO, 0, 0);

    // texcoord attribute
    glVertexArrayAttribFormat(Renderer::quadVAO, 1, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GL_FLOAT));
    glVertexArrayAttribBinding(Renderer::quadVAO, 1, 0);


    glCreateBuffers(1, &Renderer::quadVBO);
    glNamedBufferStorage(Renderer::quadVBO, array_size(quadVertices) * sizeof(float) , &quadVertices, GL_DYNAMIC_STORAGE_BIT);

}

void Renderer::ToggleCull(const Scene& scene){
    Renderer::doCull = ~Renderer::doCull;
    if (!Renderer::doCull){
        for ( const auto& [classname, entInstData] : scene.entityInstanceMap ){
            const Model& instModel = entInstData.instModel;
            glNamedBufferSubData(
                entInstData.visibleInstIndicesSSBO,
                0, 
                entInstData.instCount * sizeof(uint32_t), 
                entInstData.visibleInstIndicesBlock.data()
            );
            glNamedBufferSubData(
                instModel.drawCmdBuffer,
                0,
                instModel.drawCmdBufferVec.size() * sizeof(IndirectDrawCommand),
                instModel.drawCmdBufferVec.data()
            );
            glMemoryBarrier( GL_ALL_BARRIER_BITS );
        }

        // // set the instance counts
        // instCountSetShader->use();
        // BindDrawCmdInstCountSetBuffers(entInstData);
        // glDispatchCompute(instModel.drawCmdCount.count, 1u, 1u);

        // glMemoryBarrier( GL_ALL_BARRIER_BITS );
        // glFinish();

    }
}

void Renderer::RenderPostProcess(){
    // render
    if (mainCameraIdx == 0){
        passthroughShader->use();
    } else postProcessShader->use();
    glBindVertexArray(Renderer::quadVAO);
    glVertexArrayVertexBuffer(Renderer::quadVAO, 0, Renderer::quadVBO, 0, sizeof(float) * 4);
    glDisable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, srcFBOTexture);

    glNamedFramebufferTexture(FBO, GL_COLOR_ATTACHMENT0, dstFBOTexture, 0);
    glDrawArrays(GL_TRIANGLES, 0, 6);

}

Renderer::Renderer(float w, float h) : width(w), height(h) {
    Init(w, h);
    if (debugShader == nullptr){
        debugShader = new Shader("./resources/shader/debug.vert", "./resources/shader/debug.frag");
        debugWireShader = new Shader("./resources/shader/debug.vert", "./resources/shader/debugline.frag");

        cullShader = new ComputeShader("./resources/shader/indirect_cull_shader.comp");
        instCountSetShader = new ComputeShader("./resources/shader/indirect_instcount_set.comp");


        postProcessShader = new Shader("./resources/shader/screen.vert", "./resources/shader/screen.frag");
        passthroughShader = new Shader("./resources/shader/screen.vert", "./resources/shader/passthrough.frag");
        postProcessShader->use();
        Renderer::SetupScreenQuad();
    }

    mainCameraIdx = Renderer::allCameras.size();
    std::cout << "Main Camera Index: " << mainCameraIdx << std::endl;
    if (allCameras.size() == 0){
        Renderer::allCameras.push_back(
            Camera(w, h, glm::vec3(-3, 3.0, -3.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f, 30.0f, 0.1f, 30.0f)
        );
    } else
    {
        Renderer::allCameras.push_back(
            Camera(w, h, glm::vec3(-10.0, 20.0, -10.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f, 30.0f, 0.1f, pow( 10.0f, allCameras.size() + 1))
        );
    }
    allCameras[mainCameraIdx].setPerspectiveSize(w, h);

    Renderer::allDebugMeshes.push_back(Primitive::CreateFrustum(Renderer::allCameras[mainCameraIdx]));

    glCreateBuffers(1, &cameraMatricesUBO);
    glNamedBufferStorage(cameraMatricesUBO, 2 * sizeof(glm::mat4), NULL, GL_DYNAMIC_STORAGE_BIT);

    glCreateBuffers(1, &frustumCullDataUBO);
    glNamedBufferStorage(frustumCullDataUBO, sizeof(FrustumCullDataUBOBlock), NULL, GL_DYNAMIC_STORAGE_BIT);

    // glCreateBuffers(1, &meshPropertiesUBO);
    // glNamedBufferStorage(meshPropertiesUBO, sizeof(MeshPropertiesUBOBlock), NULL, GL_DYNAMIC_STORAGE_BIT);

    glCreateBuffers(1, &materialUniformsUBO);
    glNamedBufferStorage(materialUniformsUBO, sizeof(Material), NULL, GL_DYNAMIC_STORAGE_BIT);
    fmt::print("Material Size {}\n", sizeof(Material));

    // glCreateBuffers(1, &inDrawCmdBufferr);
    // glNamedBufferStorage(inDrawCmdBuffer,
    //                     sizeof(DrawElementsIndirectCommand) * commands.size(),
    //                     (const void *)commands.data(),
    //                     GL_DYNAMIC_STORAGE_BIT);

}

void Renderer::Init(float w, float h){
    glViewport(0, 0, w, h);
    CreateVAO();
    CreateInstanceVAO();
    CreateDebugVAO();
    CreateParticleVAO();
    CreateFBO(w, h);
    CreateRBO(w, h);

    if(glCheckNamedFramebufferStatus(FBO, GL_FRAMEBUFFER) != GL_FRAMEBUFFER_COMPLETE){
        std::cout << "ERROR::FRAMEBUFFER:: Framebuffer is not complete!" << std::endl;
        assert(false);
    } else {
        std::cout << "Framebuffer is complete!" << std::endl;
    }


}

void Renderer::CreateVAO(){
    glCreateVertexArrays(1, &VAO);

    glEnableVertexArrayAttrib(VAO, POSITION_ATTRIB_LOC);
    glEnableVertexArrayAttrib(VAO, NORMAL_ATTRIB_LOC);
    glEnableVertexArrayAttrib(VAO, TEXCOORD_ATTRIB_LOC);

    // position attribute
    glVertexArrayAttribFormat(VAO, POSITION_ATTRIB_LOC, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
    glVertexArrayAttribBinding(VAO, POSITION_ATTRIB_LOC, 0);

    // normal attribute
    glVertexArrayAttribFormat(VAO, NORMAL_ATTRIB_LOC, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));
    glVertexArrayAttribBinding(VAO, NORMAL_ATTRIB_LOC, 0);

    // texcoord attribute
    glVertexArrayAttribFormat(VAO, TEXCOORD_ATTRIB_LOC, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, texcoords));
    glVertexArrayAttribBinding(VAO, TEXCOORD_ATTRIB_LOC, 0);
}

void Renderer::CreateInstanceVAO(){
    glCreateVertexArrays(1, &instVAO);

    glEnableVertexArrayAttrib(instVAO, POSITION_ATTRIB_LOC);
    glEnableVertexArrayAttrib(instVAO, NORMAL_ATTRIB_LOC);
    glEnableVertexArrayAttrib(instVAO, TEXCOORD_ATTRIB_LOC);

    // position attribute
    glVertexArrayAttribFormat(instVAO, POSITION_ATTRIB_LOC, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
    glVertexArrayAttribBinding(instVAO, POSITION_ATTRIB_LOC, 0);

    // normal attribute
    glVertexArrayAttribFormat(instVAO, NORMAL_ATTRIB_LOC, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));
    glVertexArrayAttribBinding(instVAO, NORMAL_ATTRIB_LOC, 0);

    // texcoord attribute
    glVertexArrayAttribFormat(instVAO, TEXCOORD_ATTRIB_LOC, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, texcoords));
    glVertexArrayAttribBinding(instVAO, TEXCOORD_ATTRIB_LOC, 0);
}



void Renderer::CreateParticleVAO(){
    glCreateVertexArrays(1, &particleVAO);

    glEnableVertexArrayAttrib(particleVAO, POSITION_ATTRIB_LOC);
    
    // position attribute
    glVertexArrayAttribFormat(instVAO, POSITION_ATTRIB_LOC, 4, GL_FLOAT, GL_FALSE, 0);
    glVertexArrayAttribBinding(instVAO, POSITION_ATTRIB_LOC, 0);
}

void Renderer::CreateDebugVAO(){
    glCreateVertexArrays(1, &debugVAO);
    // glGenVertexArrays(1, &debugVAO);
    // glBindVertexArray(debugVAO);

    glEnableVertexArrayAttrib(debugVAO, 0);
    // glEnableVertexAttribArray(0);

    // position attribute
    glVertexArrayAttribFormat(debugVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
    glVertexArrayAttribBinding(debugVAO, 0, 0);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * GL_FLOAT, (void*) 0);
}

void Renderer::CreateFBO(float w, float h){
    // if (FBO != 0){
    //     glDeleteFramebuffers(1, &FBO);
    // }
    glCreateFramebuffers(1, &FBO);
    // glGenFramebuffers(1, &FBO);
    // glBindFramebuffer(GL_FRAMEBUFFER, FBO);

    // create a color attachment texture
    // if (FBOTexture != 0){
    //     glDeleteTextures(1, &FBOTexture);
    // }

    glCreateTextures(GL_TEXTURE_2D, 1, &srcFBOTexture);
    // glGenTextures(1, &srcFBOTexture);
    // glBindTexture(GL_TEXTURE_2D, srcFBOTexture);

    // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 600, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

    glTextureParameteri(srcFBOTexture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(srcFBOTexture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTextureParameteri(srcFBOTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(srcFBOTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTextureStorage2D(srcFBOTexture, 1, GL_RGBA8, w, h); // Solution

    // glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, srcFBOTexture, 0);
    // glNamedFramebufferTexture(fbo, GL_DEPTH_ATTACHMENT, depthTex, 0);


    // glBindFramebuffer(GL_FRAMEBUFFER, 0); // Reset back to default framebuffer
    // glDeleteFramebuffers(1, &FBO); // delete the existing framebuffer

    glCreateTextures(GL_TEXTURE_2D, 1, &dstFBOTexture);
    glTextureParameteri(dstFBOTexture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(dstFBOTexture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTextureParameteri(dstFBOTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(dstFBOTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTextureStorage2D(dstFBOTexture, 1, GL_RGBA8, w, h); // Solution

    // glNamedFramebufferTexture(FBO, GL_COLOR_ATTACHMENT0, dstFBOTexture, 0);

}

void Renderer::CreateRBO(float w, float h){
    // if (RBO != 0){
    //     glDeleteRenderbuffers(1, &RBO);
    // }
    glCreateRenderbuffers(1, &RBO);
    glNamedRenderbufferStorage(RBO, GL_DEPTH24_STENCIL8, w, h);
    glNamedFramebufferRenderbuffer(FBO, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);
    // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);
    // glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // glBindFramebuffer(GL_FRAMEBUFFER, 0);
    // glDeleteRenderbuffers(1, &RBO);
}

void Renderer::Resize(float w, float h){
    glViewport(0, 0, w, h);
    this->width = w;
    this->height = h;

    allCameras[mainCameraIdx].setPerspectiveSize(w, h);

    glDeleteFramebuffers(1, &FBO);
    glDeleteTextures(1, &srcFBOTexture);
    glDeleteTextures(1, &dstFBOTexture);
    glDeleteRenderbuffers(1, &RBO);

    CreateFBO(w, h);
    CreateRBO(w, h);
}


// void Renderer::BindMesh(const Primitive& mesh){
//     glVertexArrayVertexBuffer(VAO, 0, mesh.VBO, 0, sizeof(Vertex));
//     glVertexArrayElementBuffer(VAO, mesh.EBO);
//     // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
//     // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
// }

void Renderer::BindPrimitive(const Primitive& mesh){
    glVertexArrayVertexBuffer(VAO, 0, mesh.VBO, 0, sizeof(Vertex));
    glVertexArrayElementBuffer(VAO, mesh.EBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
}


void Renderer::BindDebugMesh(const Primitive& mesh){
    glVertexArrayVertexBuffer(debugVAO, 0, mesh.VBO, 0, sizeof(glm::vec3));
    glVertexArrayElementBuffer(debugVAO, mesh.EBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
}

void Renderer::BindInstanceCullingBuffers(const EntityInstanceData& entInstData){
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CULLING_WORLD_FROM_MODEL_SSBO_BINDING, entInstData.instModelMatrixBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, CULLING_VISIBLE_INSTANCES_SSBO_BINDING, entInstData.visibleInstIndicesSSBO);
    glBindBufferBase(GL_UNIFORM_BUFFER, CULLING_FRUSTUM_CULL_DATA_UBO_BINDING, frustumCullDataUBO);
}

void Renderer::BindDrawCmdInstCountSetBuffers(const EntityInstanceData& entInstData){
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTCOUNT_SET_VISIBLE_INSTANCES_SSBO_BINDING, entInstData.visibleInstIndicesSSBO);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTCOUNT_SET_DRAW_CMD_BUFFER_SSBO_BINDING, entInstData.instModel.drawCmdBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, INSTCOUNT_SET_DRAW_CMD_COUNT_SSBO_BINDING, entInstData.instModel.drawCmdCountBuffer);
}



void Renderer::BindInstanceData(const EntityInstanceData& entInstData){
    // glVertexArrayVertexBuffer(instVAO, 1, entInstData->instArrVBO, 0, sizeof(glm::mat4));

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, NODEPRIM_PROPERTIES_SSBO_BINDING, entInstData.instModel.nodePrimPropertiesBuffer);


    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PRIMITIVE_PROPERTIES_SSBO_BINDING, entInstData.instModel.primitivePropertiesBuffer);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, VISIBLE_INSTANCES_SSBO_BINDING, entInstData.visibleInstIndicesSSBO);

    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, WORLD_FROM_MODEL_SSBO_BINDING, entInstData.instModelMatrixBuffer);


    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, NODE_PROPERTIES_SSBO_BINDING, entInstData.instModel.nodePropertiesBuffer);

    // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, MODEL_FROM_MESH_SSBO_BINDING, entInstData.instModel.meshPropertiesBuffer);

    // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, MATERIAL_PROPERTIES_SSBO_BINDING, entInstData.instModel.materialPropertiesBuffer);

    // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, TEXTURES_SSBO_BINDING, entInstData.instModel.textureArray);

}

void Renderer::BindInstancePrimitive(const Primitive& primitive){
    glVertexArrayVertexBuffer(instVAO, 0, primitive.VBO, 0, sizeof(Vertex));
    glVertexArrayElementBuffer(instVAO, primitive.EBO);
}

void Renderer::BindParticleSystem(const ParticleSystem* particleSystem){
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_POSITIONS_SSBO_BINDING, particleSystem->positionsBuffer);
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_VELOCITIES_SSBO_BINDING, particleSystem->velocityBuffer);
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_FORCES_SSBO_BINDING, particleSystem->forcesBuffer);
   //glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_DATA_SSBO_BINDING, particleSystem->particleDataBuffer);
   glBindBufferBase(GL_SHADER_STORAGE_BUFFER, PARTICLE_SYSTEM_DATA_SSBO_BINDING, particleSystem->particleSystemDataBuffer);
}

void Renderer::RenderEntities(const Scene& scene){
    scene.shaders[0].use();
    // scene.shaders[0].setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
    // scene.shaders[0].setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());

    glBindVertexArray(VAO);
    for (uint32_t entityIdx = 0; entityIdx < scene.entities.size(); ++entityIdx){
        const Entity& entity = scene.entities[entityIdx];
        if (!Renderer::doCull || (Renderer::doCull && entity.model.boundingVolume.IsOnFrustum(Camera::createFrustumFromCamera(Renderer::allCameras[0]), entity.transform))){

            for ( const Node& node : entity.model.nodes ){
                const Mesh& mesh = entity.model.meshes[node.meshIndex];
                scene.shaders.at(0).setMat4("model", entity.transform.GetModelMatrix() * node.nodeTransformMatrix);
                for (const Primitive& primitive : mesh.primitives){
                    glBindTexture(GL_TEXTURE_2D, primitive.albedoTexture);
                    glBindBufferBase(GL_UNIFORM_BUFFER, MATERIAL_UBO_BINDING, materialUniformsUBO);
                    glNamedBufferSubData(materialUniformsUBO, 0, sizeof(Material), &entity.model.materials[primitive.materialUniformsIndex]);
                    // fmt::print("Alpha cutoff: {}\n", entity.model.materials[primitive.materialUniformsIndex].alphaCutoff);

                    BindPrimitive(primitive);


                    glDrawElements(GL_TRIANGLES, primitive.indexCount, GL_UNSIGNED_INT, 0);
                }
            }
        }
    }
}

void Renderer::RenderInstancedStaticModels(const Scene& scene){
    frustumCullDataUBOBlock.frustum = Camera::createFrustumFromCamera(Renderer::allCameras[0]).ToGPUFrustum();
    frustumCullDataUBOBlock.doCull = Renderer::doCull;

    for ( const auto& [classname, entInstData] : scene.entityInstanceMap ){
        const Model& instModel = entInstData.instModel;
        if (Renderer::doCull) {
            cullShader->use();
            frustumCullDataUBOBlock.boundingVolume = entInstData.instModel.boundingVolume.ToGPUSphere();
            frustumCullDataUBOBlock.instCount = entInstData.instCount;
            glNamedBufferSubData(frustumCullDataUBO, 0, sizeof(FrustumCullDataUBOBlock), &frustumCullDataUBOBlock);
            // glNamedBufferSubData(entInstData.visibleInstIndicesSSBO, 0, sizeof(uint32_t), );
            glClearNamedBufferSubData(entInstData.visibleInstIndicesSSBO, GL_R32UI, 0, sizeof(GLuint), GL_RED_INTEGER, GL_UNSIGNED_INT, NULL);
            glMemoryBarrier(GL_ALL_BARRIER_BITS);

            BindInstanceCullingBuffers(entInstData);

            glDispatchCompute(entInstData.instCount, 1u, 1u);

            glMemoryBarrier(GL_ALL_BARRIER_BITS);
            glFinish();


                // set the instance counts
                instCountSetShader->use();
                BindDrawCmdInstCountSetBuffers(entInstData);
                glDispatchCompute(instModel.drawCmdCount.count, 1u, 1u);
                glMemoryBarrier( GL_ALL_BARRIER_BITS );
                glFinish();
        }

        uint32_t visibleInstCount = 123;
        glGetNamedBufferSubData(entInstData.visibleInstIndicesSSBO, 0, sizeof(uint32_t), &visibleInstCount);
        fmt::print("Visible Instance Count: {}\n", visibleInstCount);

        scene.shaders[1].use();
        glBindBufferBase(GL_UNIFORM_BUFFER, PROJ_VIEW_UBO_BINDING, cameraMatricesUBO);
        glBindVertexArray(instVAO);

        BindInstanceData(entInstData);

        glVertexArrayVertexBuffer(instVAO, 0, instModel.VBO, 0, sizeof(Vertex));
        glVertexArrayElementBuffer(instVAO, instModel.EBO);

        glBindBuffer(GL_DRAW_INDIRECT_BUFFER, instModel.drawCmdBuffer);
        glMultiDrawElementsIndirect(
            GL_TRIANGLES,                           // Mode
            GL_UNSIGNED_INT,                        // type
            (const void*) (0 * sizeof(IndirectDrawCommand)),        // offset
            instModel.drawCmdBufferVec.size(),      // # of calls
            0                                       // stride
        );
    }
}

void Renderer::RenderParticleSystems(const Scene& scene){
    glPolygonMode( GL_FRONT_AND_BACK, GL_POINT );
    glPointSize(10.0f);

    for (const ParticleSystem& particleSystem : scene.particleSystems){
        // scene.shaders[1].use();
        particleSystem.particleComputeShader->use();
        BindParticleSystem(&particleSystem);
        glNamedBufferSubData(
            particleSystem.particleSystemDataBuffer, 
            0,
            sizeof(ParticleSystemDataBlock), 
            &particleSystem.particleSystemDataBlock
        );
        glDispatchCompute(particleSystem.particleCount, 1u, 1u);
        glMemoryBarrier( GL_ALL_BARRIER_BITS );

        glFinish();

        particleSystem.particleShader->use();
        glBindVertexArray(particleVAO);

        glVertexArrayVertexBuffer(particleVAO, 0, particleSystem.positionsBuffer, 0, sizeof(glm::vec3));
        glVertexArrayElementBuffer(particleVAO, particleSystem.EBO);

        glDrawElements(GL_POINTS, particleSystem.particleCount, GL_UNSIGNED_INT, 0);
    }
    glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );

}

void Renderer::RenderDebugVolumes(const Scene& scene){
    { // Wire Pass
        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
        debugWireShader->use();

        glBindVertexArray(debugVAO);
        for (uint32_t cameraIdx = 0; cameraIdx < Renderer::allCameras.size(); ++cameraIdx){
            if (cameraIdx == mainCameraIdx){
                continue;
            }
            debugWireShader->setMat4("model", Renderer::allCameras[cameraIdx].GetModelMatrix());
            BindDebugMesh(Renderer::allDebugMeshes[cameraIdx]);

            glDrawElements(GL_TRIANGLES, Renderer::allDebugMeshes[cameraIdx].indexCount, GL_UNSIGNED_INT, 0);
        }
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
    }


    { // Debug Camera Polygon Pass
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        debugShader->use();

        // debugShader->setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        // debugShader->setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());
        glBindVertexArray(debugVAO);
        for (uint32_t cameraIdx = 0; cameraIdx < Renderer::allCameras.size(); ++cameraIdx){
            if (cameraIdx == mainCameraIdx){
                continue;
            }
            debugShader->setMat4("model", Renderer::allCameras[cameraIdx].GetModelMatrix());
            BindDebugMesh(Renderer::allDebugMeshes[cameraIdx]);
            // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

            glDrawElements(GL_TRIANGLES, Renderer::allDebugMeshes[cameraIdx].indexCount, GL_UNSIGNED_INT, 0);

            // glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        }
    }
}

void Renderer::Render(const Scene& scene){ // really bad, we are modifying the scene state
    // render
    ZoneScoped;
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);
        glNamedFramebufferTexture(FBO, GL_COLOR_ATTACHMENT0, srcFBOTexture, 0);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // Set up
        cameraMatricesUBOBlock.projection = Renderer::allCameras[mainCameraIdx].getProjMatrix();
        cameraMatricesUBOBlock.view = Renderer::allCameras[mainCameraIdx].getViewMatrix();

        glNamedBufferSubData(cameraMatricesUBO, 0, 2 * sizeof(glm::mat4), &cameraMatricesUBOBlock);
        glBindBufferBase(GL_UNIFORM_BUFFER, PROJ_VIEW_UBO_BINDING, cameraMatricesUBO);

        RenderEntities(scene);
        RenderInstancedStaticModels(scene);
         RenderParticleSystems(scene);
        RenderDebugVolumes(scene);
        RenderPostProcess();
        // glBindVertexArray(0);
        // glUseProgram(0);
}

void Renderer::SwitchCamera(const Scene& scene){
    mainCameraIdx = (mainCameraIdx + 1) % Renderer::allCameras.size();


    std::cout << "Now using Camera #" << mainCameraIdx << std::endl;
}
