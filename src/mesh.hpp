#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <vector>
#include <string>

#define ARRAY_COUNT(x) (sizeof(x) / sizeof((x)[0]))

// https://www.wikihow.com/Linux-How-to-Mount-Drive#:~:text=To%20mount%20a%20drive%20on,to%20mount%20and%20unmount%20drives.

class Camera;

struct Vertex{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 texcoords;
};


struct Mesh {
    // public:
    //     // Create a renderer class and compare times

        GLuint VBO, VAO, EBO;
        static void GenerateBuffers(Mesh* mesh, const std::vector<Vertex> &vertices, const std::vector<uint>& indices);
        static void GenerateDebugBuffers(Mesh* mesh, const std::vector<glm::vec3> &vertices, const std::vector<uint>& indices);

        GLsizei indexCount;

        static Mesh CreateCube();
        static Mesh CreateFrustum(const Camera& camera);
        // static Mesh CreatePyramid();

    //     std::vector<float> vertices;
    //     std::vector<unsigned int> indices;
    //     // std::vector<float> tricolors;
    //     virtual void calculateVertices();
    //     void initializeBuffers();

    //     // virtual void init();

    // public:
    //     std::string type;
    //     glm::vec4 color;
    //     glm::mat4 modelMat;

    // // private:
    // //     void drawTriangle3D();
    // //     void pushTriangle3D(float verts[]);


    //     Mesh();
    //     ~Mesh();

    //     Mesh(const Mesh& other) = delete;
    //     Mesh& operator=(const Mesh& other) = delete;

    //     Mesh(Mesh&& other);
    //     Mesh& operator=(Mesh&& other);

    //     void deleteBuffers();


    //     void translate(glm::vec3 delta){
    //         modelMat = glm::translate(modelMat, delta);
    //     };
    //     void rotate(float angle, glm::vec3 axis){
    //         modelMat = glm::rotate(modelMat, angle, axis);
    //     };
    //     void scale(glm::vec3 factor){
    //         modelMat = glm::scale(modelMat, factor);
    //     };

    //     void reset() {
    //         modelMat = glm::mat4(1.0);

    //     };

    //     void printVertices();

    //     virtual void render(); // I believe we can optimize this from asg3
};
#endif