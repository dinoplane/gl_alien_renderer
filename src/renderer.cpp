#include <renderer.hpp>
#include <camera.hpp>
#include <scene.hpp>
#include <shader_s.hpp>
#include <mesh.hpp>
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


Renderer::Renderer(float w, float h) : width(w), height(h) {
    Init();
    if (debugShader == nullptr){
        debugShader = new Shader("./resources/shader/debug.vert", "./resources/shader/debug.frag");
        debugWireShader = new Shader("./resources/shader/debug.vert", "./resources/shader/debugline.frag");
    }

    mainCameraIdx = Renderer::allCameras.size();
    std::cout << "Main Camera Index: " << mainCameraIdx << std::endl;
    if (allCameras.size() == 0){
        Renderer::allCameras.push_back(
            Camera(w, h, glm::vec3(0.51, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0), 0.0f, 0.0f, 30.0f, 0.1f, 10.0f)
        );
    } else
    {
        Renderer::allCameras.push_back(
            Camera(w, h, glm::vec3(0.0, 20.0, 0.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f, 30.0f, 0.1f, pow( 10.0f, allCameras.size() + 1))
        );
    }

    Renderer::allDebugMeshes.push_back(Mesh::CreateFrustum(Renderer::allCameras[mainCameraIdx]));

}

void Renderer::Init(){

    CreateVAO();

    CreateDebugVAO();

    CreateFBO();

    CreateRBO();

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

void Renderer::CreateFBO(){
    glCreateFramebuffers(1, &FBO);
    // glGenFramebuffers(1, &FBO);
    // glBindFramebuffer(GL_FRAMEBUFFER, FBO);

    // create a color attachment texture

    glCreateTextures(GL_TEXTURE_2D, 1, &FBOTexture);
    // glGenTextures(1, &FBOTexture);
    // glBindTexture(GL_TEXTURE_2D, FBOTexture);

    // glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, 800, 600, 0, GL_RGB, GL_UNSIGNED_BYTE, NULL);

    glTextureParameteri(FBOTexture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTextureParameteri(FBOTexture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

    glTextureStorage2D(FBOTexture, 1, GL_RGBA8, width, height);
    // glTextureSubImage2D(FBOTexture, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);
    // glTexStorage2D(GL_TEXTURE_2D, 1, GL_RGBA8, width, height);
    // glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, pixels);

    glNamedFramebufferTexture(FBO, GL_COLOR_ATTACHMENT0, FBOTexture, 0);
    // glFramebufferTexture2D(GL_FRAMEBUFFER, GL_COLOR_ATTACHMENT0, GL_TEXTURE_2D, FBOTexture, 0);
    // glNamedFramebufferTexture(fbo, GL_DEPTH_ATTACHMENT, depthTex, 0);




    // glBindFramebuffer(GL_FRAMEBUFFER, 0); // Reset back to default framebuffer
    // glDeleteFramebuffers(1, &FBO); // delete the existing framebuffer

}

void Renderer::CreateRBO(){
    glCreateRenderbuffers(1, &RBO);
    // glGenRenderbuffers(1, &RBO);
    // glBindRenderbuffer(GL_RENDERBUFFER, RBO);

    glNamedRenderbufferStorage(RBO, GL_DEPTH24_STENCIL8, width, height);
    // glRenderbufferStorage(GL_RENDERBUFFER, GL_DEPTH24_STENCIL8, 800, 600);

    glNamedFramebufferRenderbuffer(FBO, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);
    // glFramebufferRenderbuffer(GL_FRAMEBUFFER, GL_DEPTH_STENCIL_ATTACHMENT, GL_RENDERBUFFER, RBO);
    // glBindRenderbuffer(GL_RENDERBUFFER, 0);
    // glBindFramebuffer(GL_FRAMEBUFFER, 0);
    // glDeleteRenderbuffers(1, &RBO);
}

void Renderer::Resize(int w, int h){
    glViewport(0, 0, w, h);
    this->width = w;
    this->height = h;

    allCameras[mainCameraIdx].setPerspectiveSize(w, h);
}


void Renderer::BindMesh(Mesh * mesh){
    glVertexArrayVertexBuffer(VAO, 0, mesh->VBO, 0, sizeof(Vertex));
    glVertexArrayElementBuffer(VAO, mesh->EBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
}

void Renderer::BindDebugMesh(Mesh * mesh){
    glVertexArrayVertexBuffer(debugVAO, 0, mesh->VBO, 0, sizeof(glm::vec3));
    glVertexArrayElementBuffer(debugVAO, mesh->EBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
}

void Renderer::Render(Scene* scene){ // really bad, we are modifying the scene state
    // render
        glBindFramebuffer(GL_FRAMEBUFFER, FBO);

        glEnable(GL_DEPTH_TEST);
        glDepthFunc(GL_LESS);
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        scene->shaders[0].use();
        scene->shaders[0].setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        scene->shaders[0].setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());
        glBindVertexArray(VAO);
        for (uint entityIdx = 0; entityIdx < scene->entities.size(); ++entityIdx){
            scene->shaders.at(0).setMat4("model", scene->entities[entityIdx].transform.GetModelMatrix());
            BindMesh(&scene->entities[entityIdx].mesh);
            glDrawElements(GL_TRIANGLES, scene->entities[entityIdx].mesh.indexCount, GL_UNSIGNED_INT, 0);
        }


        glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
        debugWireShader->use();
        debugWireShader->setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        debugWireShader->setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());
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
        debugShader->setMat4("projection", Renderer::allCameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        debugShader->setMat4("view", Renderer::allCameras[mainCameraIdx].getViewMatrix());
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
        // glBindVertexArray(0);
        // glUseProgram(0);
}

void Renderer::SwitchCamera(const Scene& scene){
    mainCameraIdx = (mainCameraIdx + 1) % Renderer::allCameras.size();


    std::cout << "Now using Camera #" << mainCameraIdx << std::endl;
}
