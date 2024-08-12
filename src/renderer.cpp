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


Renderer::Renderer(){
    Init();
    if (debugShader == nullptr){
        debugShader = new Shader("./resources/shader/debug.vert", "./resources/shader/debug.frag");
        debugWireShader = new Shader("./resources/shader/debug.vert", "./resources/shader/debugline.frag");
    }

    mainCameraIdx = Renderer::allCameras.size();
    std::cout << "Main Camera Index: " << mainCameraIdx << std::endl;
    if (allCameras.size() == 0){
        Renderer::allCameras.push_back(
            Camera(800.0, 600.0, glm::vec3(0.51, 0.0, 0.0), glm::vec3(0.0, 1.0, 0.0), 0.0f, 0.0f, 30.0f, 0.1f, 10.0f)
        );
    } else
    {
        Renderer::allCameras.push_back(
            Camera(800.0, 600.0, glm::vec3(0.0, 20.0, 0.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f, 30.0f, 0.1f, pow( 10.0f, allCameras.size() + 1))
        );
    }

    Renderer::allDebugMeshes.push_back(Mesh::CreateFrustum(Renderer::allCameras[mainCameraIdx]));

}

void Renderer::Init(){
    // glGenVertexArrays(1, &VAO);
    BindVAO();

    // glGenVertexArrays(1, &debugVAO);
    BindDebugVAO();
}

void Renderer::BindVAO(){
    glCreateVertexArrays(1, &VAO);
    // glBindVertexArray(VAO);

    glEnableVertexArrayAttrib(VAO, 0);
    glEnableVertexArrayAttrib(VAO, 1);
    glEnableVertexArrayAttrib(VAO, 2);

    // position attribute
    glVertexArrayAttribFormat(VAO, 0, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
    glVertexArrayAttribBinding(VAO, 0, 0);
    // glEnableVertexAttribArray(0);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, position));

    // std::cout << "Position offset: " << offsetof(Vertex, position) << std::endl;

    // normal attribute
    glVertexArrayAttribFormat(VAO, 1, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));
    glVertexArrayAttribBinding(VAO, 1, 0);
    // glEnableVertexAttribArray(1);
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, normal));

    // std::cout << "Normal offset: " << offsetof(Vertex, normal) << std::endl;

    // texcoord attribute
    glVertexArrayAttribFormat(VAO, 2, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, texcoords));
    glVertexArrayAttribBinding(VAO, 2, 0);
    // glEnableVertexAttribArray(2);
    // glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, sizeof(Vertex), (void*) offsetof(Vertex, texcoords));

    // std::cout << "Texcoords offset: " << offsetof(Vertex, texcoords) << std::endl;

}

void Renderer::BindDebugVAO(){
    glCreateVertexArrays(1, &debugVAO);
    // glBindVertexArray(debugVAO);

    glEnableVertexArrayAttrib(debugVAO, 0);

    // position attribute
    glVertexArrayAttribFormat(debugVAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
    glVertexArrayAttribBinding(debugVAO, 0, 0);
    // glEnableVertexAttribArray(0);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 3 * GL_FLOAT, (void*) 0);
}


void Renderer::Resize(int width, int height){
    glViewport(0, 0, width, height);
    allCameras[mainCameraIdx].setPerspectiveSize(width, height);
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
