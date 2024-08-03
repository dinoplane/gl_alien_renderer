#include <renderer.hpp>
#include <camera.hpp>
#include <scene.hpp>
#include <shader_s.hpp>
#include <mesh.hpp>

#include <glad/glad.h>

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

        scene->shaders[0].use();
        scene->shaders[0].setMat4("projection", scene->cameras[0].getProjMatrix()); // TODO : Profile this
        scene->shaders[0].setMat4("view", scene->cameras[0].getViewMatrix());
        // scene->shaders[0].setMat4("projection", glm::identity<glm::mat4>()); // TODO : Profile this
        // scene->shaders[0].setMat4("view", glm::identity<glm::mat4>());
        scene->shaders[0].setMat4("model", glm::identity<glm::mat4>());
        glBindVertexArray(scene->meshes[0].VAO);
        glDrawElements(GL_TRIANGLES, scene->meshes[0].indexCount, GL_UNSIGNED_INT, 0);
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


