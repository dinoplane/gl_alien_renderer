#include <scenedata.hpp>

#include <util.h>

#include <shader_s.hpp>
#include <mesh.hpp>
#include <transform.hpp>
#include <entity.hpp>
#include <camera.hpp>

#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <glm/gtx/string_cast.hpp>

MeshData MeshData::CreateCube(){ // Need to change this later
    return {
                {
                    // counterclockwise
                    // positions                   // normals                     // texcoords

                    // top face
                    {{-0.5f, -0.5f,  0.5f},        { 0.0f,  0.0f,  1.0f},         {0.0f, 0.0f}},
                    {{ 0.5f, -0.5f,  0.5f},        { 0.0f,  0.0f,  1.0f},         {1.0f, 0.0f}},
                    {{ 0.5f,  0.5f,  0.5f},        { 0.0f,  0.0f,  1.0f},         {1.0f, 1.0f}},
                    {{-0.5f,  0.5f,  0.5f},        { 0.0f,  0.0f,  1.0f},         {0.0f, 1.0f}},

                    // right face
                    {{ 0.5f, -0.5f, -0.5f},        { 1.0f,  0.0f,  0.0f},         {0.0f, 0.0f}},
                    {{ 0.5f,  0.5f, -0.5f},        { 1.0f,  0.0f,  0.0f},         {1.0f, 0.0f}},
                    {{ 0.5f,  0.5f,  0.5f},        { 1.0f,  0.0f,  0.0f},         {1.0f, 1.0f}},
                    {{ 0.5f, -0.5f,  0.5f},        { 1.0f,  0.0f,  0.0f},         {0.0f, 1.0f}},

                    // bot face
                    {{-0.5f,  0.5f, -0.5f},        { 0.0f,  0.0f, -1.0f},         {0.0f, 0.0f}},
                    {{ 0.5f,  0.5f, -0.5f},        { 0.0f,  0.0f, -1.0f},         {1.0f, 0.0f}},
                    {{ 0.5f, -0.5f, -0.5f},        { 0.0f,  0.0f, -1.0f},         {1.0f, 1.0f}},
                    {{-0.5f, -0.5f, -0.5f},        { 0.0f,  0.0f, -1.0f},         {0.0f, 1.0f}},

                    // front face
                    {{ 0.5f,  0.5f, -0.5f},        { 0.0f,  1.0f,  0.0f},         {0.0f, 0.0f}},
                    {{-0.5f,  0.5f, -0.5f},        { 0.0f,  1.0f,  0.0f},         {1.0f, 0.0f}},
                    {{-0.5f,  0.5f,  0.5f},        { 0.0f,  1.0f,  0.0f},         {1.0f, 1.0f}},
                    {{ 0.5f,  0.5f,  0.5f},        { 0.0f,  1.0f,  0.0f},         {0.0f, 1.0f}},

                    // left face
                    {{-0.5f,  0.5f, -0.5f},        {-1.0f,  0.0f,  0.0f},         {0.0f, 0.0f}},
                    {{-0.5f, -0.5f, -0.5f},        {-1.0f,  0.0f,  0.0f},         {1.0f, 0.0f}},
                    {{-0.5f, -0.5f,  0.5f},        {-1.0f,  0.0f,  0.0f},         {1.0f, 1.0f}},
                    {{-0.5f,  0.5f,  0.5f},        {-1.0f,  0.0f,  0.0f},         {0.0f, 1.0f}},

                    // back face
                    {{-0.5f, -0.5f, -0.5f},        { 0.0f, -1.0f,  0.0f},         {0.0f, 0.0f}},
                    {{ 0.5f, -0.5f, -0.5f},        { 0.0f, -1.0f,  0.0f},         {1.0f, 0.0f}},
                    {{ 0.5f, -0.5f,  0.5f},        { 0.0f, -1.0f,  0.0f},         {1.0f, 1.0f}},
                    {{-0.5f, -0.5f,  0.5f},        { 0.0f, -1.0f,  0.0f},         {0.0f, 1.0f}},
                },
                {  // note that we start from 0!
                    // top
                        0, 1, 2, 2, 3, 0,
                    // right
                        4, 5, 6, 6, 7, 4,
                    // bot
                        8, 9, 10, 10, 11, 8,
                    // front
                        12, 13, 14, 14, 15, 12,
                    // left
                        16, 17, 18, 18, 19, 16,
                    // back
                        20, 21, 22, 22, 23, 20
                }
            };
}

SceneData SceneData::GenerateDefaultScene(){
    SceneData retSceneData;

    Transform transform;

    for (int i = 0; i < 5; ++i){
        for (int j = 0; j < 5; ++j){
            for (int k = 0; k < 5; ++k){
                transform.SetPosition(glm::vec3(i * 2.0, j * 2.0, k * 2.0));
                transform.SetRotation(glm::vec3(0.0, 0.0, 0.0));
                transform.SetScale(glm::vec3(1.0, 1.0, 1.0));
                retSceneData.entitiesData.push_back({MeshData::CreateCube(), transform}); // TODO Im pushing back a lot of data...
            }
        }
    }

    // transform.SetPosition(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetRotation(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetScale(glm::vec3(1.0, 1.0, 10.0));
    // // transform.GetModelMatrix();
    // retSceneData.entities.push_back({Mesh::CreateCube(), transform});


    // transform.SetPosition(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetRotation(glm::vec3(0.0, 90.0, 0.0));
    // transform.SetScale(glm::vec3(1.0, 10.0, 1.0));
    // retSceneData.entities.push_back({Mesh::CreateCube(), transform});

    retSceneData.shaderPaths.push_back({"./resources/shader/base.vert", "./resources/shader/base.frag"});
    retSceneData.shaderPaths.push_back({"./resources/shader/debug.vert", "./resources/shader/debug.frag"});


    retSceneData.cameraData.push_back(Camera(800.0, 600.0, glm::vec3(0.0, 20.0, 0.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f));
    retSceneData.cameraData.push_back(Camera(800.0, 600.0));
    return retSceneData;
}