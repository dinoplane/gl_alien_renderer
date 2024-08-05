#include <renderer.hpp>
#include <camera.hpp>
#include <scene.hpp>
#include <shader_s.hpp>
#include <mesh.hpp>
#include <entity.hpp>

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <iostream>

void Renderer::Init(){

}

void Renderer::Resize(const Camera& camera){
// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
// void framebuffer_size_callback(GLFWwindow* window, int width, int height)
// {
//     // make sure the viewport matches the new window dimensions; note that width and
//     // height will be significantly larger than specified on retina displays.
//     glViewport(0, 0, width, height);
//     camera->setPerspectiveSize(width, height);
// }

}

void Renderer::Render(Scene* scene){ // really bad, we are modifying the scene state
    // render

        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // scene->shaders[0].use();
        // scene->shaders[0].setMat4("projection", scene->cameras[0].getProjMatrix()); // TODO : Profile this
        // scene->shaders[0].setMat4("view", scene->cameras[0].getViewMatrix());
        // scene->shaders[0].setMat4("model", scene->entities[0].transform.GetModelMatrix());

        // // std::cout << "Model Matrix : " << glm::to_string(scene->entities[0].transform.GetModelMatrix()) << std::endl;
        // // scene->shaders[0].setMat4("projection", glm::identity<glm::mat4>()); // TODO : Profile this
        // // scene->shaders[0].setMat4("view", glm::identity<glm::mat4>());
        // // scene->shaders[0].setMat4("model", glm::identity<glm::mat4>());
        // glBindVertexArray(scene->entities[0].mesh.VAO);
        // glDrawElements(GL_TRIANGLES, scene->entities[0].mesh.indexCount, GL_UNSIGNED_INT, 0);

        scene->shaders[0].use();
        scene->shaders[0].setMat4("projection", scene->cameras[mainCameraIdx].getProjMatrix()); // TODO : Profile this
        scene->shaders[0].setMat4("view", scene->cameras[mainCameraIdx].getViewMatrix());
        // std::cout << "Shader Size: " << scene->shaders.size() << std::endl;

        for (uint entityIdx = 0; entityIdx < scene->entities.size(); ++entityIdx){
            scene->shaders.at(0).setMat4("model", scene->entities[entityIdx].transform.GetModelMatrix());
            // std::cout << "EntityIndex " << entityIdx << std::endl;

            // std::cout << "Model Matrix : " << glm::to_string(scene->entities[entityIdx].transform.GetModelMatrix()) << std::endl;
            glBindVertexArray(scene->entities[entityIdx].mesh.VAO);
            glDrawElements(GL_TRIANGLES, scene->entities[entityIdx].mesh.indexCount, GL_UNSIGNED_INT, 0);
        }
        // Another gl draw elements for the camera debug

        // ourShader.use();
        // ourShader.setMat4("projection", projection);
        // ourShader.setMat4("view", view);

        // ourShader.setInt("selected", 1);
        // ourShader.setMat4("model", floor.modelMat);
        // // seaShader.setInt("uIsFloor", 0);

        // seafloorShader.use();
        // seafloorShader.setMat4("projection", projection);
        // seafloorShader.setMat4("view", view);
        // seafloorShader.setInt("selected", 0);
        // seafloorShader.setFloat("uTime", glfwGetTime());
        // seafloorShader.setMat4("model", floor.modelMat);
        // seafloorShader.setVec3("uLightPos", lightPos);
        // seafloorShader.setVec3("uViewPos", camera.position);
        // floor.render();

        // ourShader.use();
        // ourShader.setMat4("model", xcoord.modelMat);
        // xcoord.render();

        // ourShader.use();
        // ourShader.setMat4("model", ycoord.modelMat);
        // ycoord.render();

        // ourShader.use();
        // ourShader.setMat4("model", zcoord.modelMat);
        // zcoord.render();

        // ourShader.use();
        // ourShader.setMat4("model", light.modelMat);
        // light.render();

        // seaShader.use();
        // seaShader.setMat4("projection", projection);
        // seaShader.setMat4("view", view);
        // seaShader.setInt("selected", 0);
        // seaShader.setFloat("uTime", glfwGetTime());
        // seaShader.setMat4("model", sea.modelMat);
        // seaShader.setVec3("uLightPos", lightPos);
        // seaShader.setVec3("uViewPos", camera.position);
        // sea.render();


}

void Renderer::SwitchCamera(const Scene& scene){
    mainCameraIdx = (mainCameraIdx + 1) % scene.cameras.size();


    std::cout << "Now using Camera #" << mainCameraIdx << std::endl;
}
