#include <renderer.hpp>
#include <camera.hpp>
#include <scene.hpp>
#include <shader_s.hpp>
#include <mesh.hpp>
#include <shader_c.hpp>
#include <entity.hpp>
#include <util.h>


#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>

std::vector<Camera> Renderer::allCameras;
std::vector<Mesh> Renderer::allDebugMeshes;
Shader* Renderer::debugShader = nullptr;
Shader* Renderer::debugWireShader = nullptr;
Shader* Renderer::postProcessShader = nullptr;
Shader* Renderer::passthroughShader = nullptr;
ComputeShader* Renderer::cullShader = nullptr;
GLuint Renderer::quadVAO = 0;
GLuint Renderer::quadVBO = 0;
uint Renderer::doCull = 0;


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

        cullShader = new ComputeShader("./resources/shader/cull_shader.comp");


        postProcessShader = new Shader("./resources/shader/screen.vert", "./resources/shader/screen.frag");
        passthroughShader = new Shader("./resources/shader/screen.vert", "./resources/shader/passthrough.frag");
        postProcessShader->use();
        Renderer::SetupScreenQuad();
    }

    mainCameraIdx = Renderer::allCameras.size();
    std::cout << "Main Camera Index: " << mainCameraIdx << std::endl;
    if (allCameras.size() == 0){
        Renderer::allCameras.push_back(
            Camera(w, h, glm::vec3(-3, 3.0, -3.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f, 30.0f, 0.1f, 10.0f)
        );
    } else
    {
        Renderer::allCameras.push_back(
            Camera(w, h, glm::vec3(-10.0, 20.0, -10.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f, 30.0f, 0.1f, pow( 10.0f, allCameras.size() + 1))
        );
    }
    allCameras[mainCameraIdx].setPerspectiveSize(w, h);

    Renderer::allDebugMeshes.push_back(Mesh::CreateFrustum(Renderer::allCameras[mainCameraIdx]));

    glCreateBuffers(1, &cameraMatricesUBO);
    glNamedBufferStorage(cameraMatricesUBO, 2 * sizeof(glm::mat4), NULL, GL_DYNAMIC_STORAGE_BIT);

    glCreateBuffers(1, &frustumCullDataUBO);
    glNamedBufferStorage(frustumCullDataUBO, sizeof(FrustumCullDataUBOBlock), NULL, GL_DYNAMIC_STORAGE_BIT);

}

void Renderer::Init(float w, float h){
    glViewport(0, 0, w, h);
    CreateVAO();
    CreateInstanceVAO();
    CreateDebugVAO();
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
    // glGenVertexArrays(1, &VAO);
    // glBindVertexArray(VAO);

    glEnableVertexArrayAttrib(VAO, 0);
    glEnableVertexArrayAttrib(VAO, 1);
    glEnableVertexArrayAttrib(VAO, 2);
    // glEnableVertexAttribArray(0);
    // glEnableVertexAttribArray(1);
    // glEnableVertexAttribArray(2);

    // position attribute
    glVertexArrayAttribFormat(VAO, 0, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
    glVertexArrayAttribBinding(VAO, 0, 0);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, position));

    // normal attribute
    glVertexArrayAttribFormat(VAO, 1, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));
    glVertexArrayAttribBinding(VAO, 1, 0);
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, normal));

    // texcoord attribute
    glVertexArrayAttribFormat(VAO, 2, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, texcoords));
    glVertexArrayAttribBinding(VAO, 2, 0);
    // glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, texcoords));
}

void Renderer::CreateInstanceVAO(){
    glCreateVertexArrays(1, &instVAO);

    glEnableVertexArrayAttrib(instVAO, 0);
    glEnableVertexArrayAttrib(instVAO, 1);
    glEnableVertexArrayAttrib(instVAO, 2);

    // glEnableVertexArrayAttrib(instVAO, 3);
    // glEnableVertexArrayAttrib(instVAO, 4);
    // glEnableVertexArrayAttrib(instVAO, 5);
    // glEnableVertexArrayAttrib(instVAO, 6);


    // position attribute
    glVertexArrayAttribFormat(instVAO, 0, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
    glVertexArrayAttribBinding(instVAO, 0, 0);

    // normal attribute
    glVertexArrayAttribFormat(instVAO, 1, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));
    glVertexArrayAttribBinding(instVAO, 1, 0);

    // texcoord attribute
    glVertexArrayAttribFormat(instVAO, 2, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, texcoords));
    glVertexArrayAttribBinding(instVAO, 2, 0);

    // modelmatrix attribute
    // glVertexArrayAttribFormat(instVAO, 3, 4, GL_FLOAT, GL_FALSE, 0 * sizeof(glm::vec4));
    // glVertexArrayAttribBinding(instVAO, 3, 1);
    // glVertexArrayAttribFormat(instVAO, 4, 4, GL_FLOAT, GL_FALSE, 1 * sizeof(glm::vec4));
    // glVertexArrayAttribBinding(instVAO, 4, 1);
    // glVertexArrayAttribFormat(instVAO, 5, 4, GL_FLOAT, GL_FALSE, 2 * sizeof(glm::vec4));
    // glVertexArrayAttribBinding(instVAO, 5, 1);
    // glVertexArrayAttribFormat(instVAO, 6, 4, GL_FLOAT, GL_FALSE, 3 * sizeof(glm::vec4));
    // glVertexArrayAttribBinding(instVAO, 6, 1);

    // glVertexArrayBindingDivisor(instVAO, 1, 1);
    // glVertexArrayBindingDivisor(instVAO, 4, 1);
    // glVertexArrayBindingDivisor(instVAO, 5, 1);
    // glVertexArrayBindingDivisor(instVAO, 6, 1);

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


void Renderer::BindMesh(const Mesh * mesh){
    glVertexArrayVertexBuffer(VAO, 0, mesh->VBO, 0, sizeof(Vertex));
    glVertexArrayElementBuffer(VAO, mesh->EBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
}

void Renderer::BindDebugMesh(const Mesh * mesh){
    glVertexArrayVertexBuffer(debugVAO, 0, mesh->VBO, 0, sizeof(glm::vec3));
    glVertexArrayElementBuffer(debugVAO, mesh->EBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
}

void Renderer::BindInstanceCullingBuffers(const EntityInstanceData* entInstData){
    // Input --------------------------------------
    // Bind the bounding volumes
    // glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, entInstData->instBoundingVolumeBuffer);

    // bind the model matrices
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, entInstData->instModelMatrixBuffer);


    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, entInstData->instMeshRenderedBuffer);


    // cullShader->use()


    // Input --------------------------------------


}


void Renderer::BindInstanceMesh(const EntityInstanceData* entInstData){
    glVertexArrayVertexBuffer(instVAO, 0, entInstData->instMesh.VBO, 0, sizeof(Vertex));
    glVertexArrayElementBuffer(instVAO, entInstData->instMesh.EBO);

    // glVertexArrayVertexBuffer(instVAO, 1, entInstData->instArrVBO, 0, sizeof(glm::mat4));
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 4, entInstData->instModelMatrixBuffer);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 5, entInstData->instMeshRenderedBuffer);

}

void Renderer::Render(const Scene* scene){ // really bad, we are modifying the scene state
    // render
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
        glBindBufferBase(GL_UNIFORM_BUFFER, 3, cameraMatricesUBO);

        scene->shaders[0].use();
        // scene->shaders[0].setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        // scene->shaders[0].setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());

        glBindVertexArray(VAO);
        for (uint entityIdx = 0; entityIdx < scene->entities.size(); ++entityIdx){
            if (scene->entities[entityIdx].mesh.boundingVolume->IsOnFrustum(Camera::createFrustumFromCamera(Renderer::allCameras[0]), scene->entities[entityIdx].transform)){
                scene->shaders.at(0).setMat4("model", scene->entities[entityIdx].transform.GetModelMatrix());
                BindMesh(&scene->entities[entityIdx].mesh);
                glDrawElements(GL_TRIANGLES, scene->entities[entityIdx].mesh.indexCount, GL_UNSIGNED_INT, 0);
            }
        }

        // cullShader->

        // glBindVertexArray(instVAO);
        // for (uint entityIdx = 0; entityIdx < scene->entityInstanceMap.size(); ++entityIdx){
        //     // scene->shaders.at(0).setMat4("model", scene->entityInstanceMap[entityIdx].transform.GetModelMatrix());
        //     // What I'm about to do justifies the need to separate static state from dynamic state

        //     BindInstanceMesh(&scene->entityInstanceMap[entityIdx]);
        //     glDrawElementsInstanced(GL_TRIANGLES, scene->entityInstanceMap[entityIdx].instMesh.indexCount, GL_UNSIGNED_INT, 0, scene->entityInstanceMap[entityIdx].instCount);
        // }


        // scene->shaders[1].setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        // scene->shaders[1].setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());

        frustumCullDataUBOBlock.frustum = Camera::createFrustumFromCamera(Renderer::allCameras[0]).ToGPUFrustum();
        frustumCullDataUBOBlock.doCull = doCull;
        
        for ( const auto& [classname, entInstData] : scene->entityInstanceMap ){
            // scene->shaders.at(0).setMat4("model", scene->entityInstanceMap[entityIdx].transform.GetModelMatrix());
            // What I'm about to do justifies the need to separate static state from dynamic state
            frustumCullDataUBOBlock.boundingVolume = entInstData.instMesh.boundingVolume->ToGPUSphere();
            frustumCullDataUBOBlock.instCount = entInstData.instCount;
            //fmt::print("Instance Count: {}\n", entInstData.instCount);
            glNamedBufferSubData(frustumCullDataUBO, 0, sizeof(FrustumCullDataUBOBlock), &frustumCullDataUBOBlock);
            glBindBufferBase(GL_UNIFORM_BUFFER, 6, frustumCullDataUBO);

            cullShader->use();
            BindInstanceCullingBuffers(&entInstData);
            glDispatchCompute(entInstData.instCount, 1u, 1u);

            glMemoryBarrier( GL_ALL_BARRIER_BITS );
            glFinish();


            scene->shaders[1].use();
            glBindBufferBase(GL_UNIFORM_BUFFER, 3, cameraMatricesUBO);

            glBindVertexArray(instVAO);


            BindInstanceMesh(&entInstData);
            glDrawElementsInstanced(GL_TRIANGLES, entInstData.instMesh.indexCount, GL_UNSIGNED_INT, 0, entInstData.instCount);
        }


        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
        debugWireShader->use();

        // debugWireShader->setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        // debugWireShader->setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());
        glBindVertexArray(debugVAO);
        for (uint cameraIdx = 0; cameraIdx < Renderer::allCameras.size(); ++cameraIdx){
            if (cameraIdx == mainCameraIdx){
                continue;
            }
            debugWireShader->setMat4("model", Renderer::allCameras[cameraIdx].GetModelMatrix());
            BindDebugMesh(&Renderer::allDebugMeshes[cameraIdx]);

            glDrawElements(GL_TRIANGLES, Renderer::allDebugMeshes[cameraIdx].indexCount, GL_UNSIGNED_INT, 0);
        }
        glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );


        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
        debugShader->use();

        // debugShader->setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        // debugShader->setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());
        glBindVertexArray(debugVAO);
        for (uint cameraIdx = 0; cameraIdx < Renderer::allCameras.size(); ++cameraIdx){
            if (cameraIdx == mainCameraIdx){
                continue;
            }
            debugShader->setMat4("model", Renderer::allCameras[cameraIdx].GetModelMatrix());
            BindDebugMesh(&Renderer::allDebugMeshes[cameraIdx]);
            // glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );

            glDrawElements(GL_TRIANGLES, Renderer::allDebugMeshes[cameraIdx].indexCount, GL_UNSIGNED_INT, 0);

            // glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
        }

        RenderPostProcess();
        // glBindVertexArray(0);
        // glUseProgram(0);
}

void Renderer::SwitchCamera(const Scene& scene){
    mainCameraIdx = (mainCameraIdx + 1) % Renderer::allCameras.size();


    std::cout << "Now using Camera #" << mainCameraIdx << std::endl;
}
