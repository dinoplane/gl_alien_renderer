#ifndef MESH_H
#define MESH_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>

#include <util.h>
#include <volume.hpp>
#include <material.hpp>

#include <vector>
#include <string>



// https://www.wikihow.com/Linux-How-to-Mount-Drive#:~:text=To%20mount%20a%20drive%20on,to%20mount%20and%20unmount%20drives.

class Camera;
class Renderer;

struct Primitive;

struct IndirectDrawCommand {
	uint32_t count;
	uint32_t instanceCount;
	uint32_t firstIndex;
    uint32_t baseVertex;
	uint32_t baseInstance;
};


struct Mesh {
    GLuint drawsBuffer;
    std::vector<Primitive> primitives;
    uint32_t startPrimIdx;
    uint32_t primCount;
    Sphere boundingVolume;

    // Sphere* boundingVolume;
    static Mesh CreateCube();
    static Mesh CreateFrustum(const Camera& camera);
    static Mesh CreatePyramid();
};

struct Primitive {
    // public:
    //     // Create a renderer class and compare times
        IndirectDrawCommand draw;
        GLuint VBO, VAO, EBO; // keep the VAO in there, just in case.
        
        uint32_t vertexStartIdx;
        uint32_t vertexCount;
        uint32_t indexStartIdx;
        uint32_t indexCount;

        size_t materialUniformsIndex;
        GLuint albedoTexture;

        static void GenerateBuffers(Primitive* primitive, const std::vector<Vertex> &vertices, const std::vector<uint32_t>& indices);
        // static void Rebind(Mesh*primitive, Renderer* renderer);
        static void GenerateDebugBuffers(Primitive* primitive, const std::vector<glm::vec3> &vertices, const std::vector<uint32_t>& indices);
        // static void RebindDebug(Mesh*primitive);


        static Primitive CreateCube();
        static Primitive CreateFrustum(const Camera& camera);
        static Primitive CreatePyramid();
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