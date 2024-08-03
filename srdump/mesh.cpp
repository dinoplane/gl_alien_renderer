#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glm/gtx/string_cast.hpp>

#include<mesh.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>


Mesh::Mesh(){


    modelMat = glm::mat4(1.0f);
    calculateVertices();
    initializeBuffers();

}

Mesh::~Mesh(){
    deleteBuffers();
}

Mesh::Mesh(Mesh&& other) : VBO(other.VBO), VAO(other.VAO), EBO(other.EBO), vertices(std::move(other.vertices)), indices(std::move(other.indices)){

    other.VBO = other.VAO = other.EBO = 0;
}

Mesh& Mesh::operator=(Mesh&& other){
    if (this != &other){
        deleteBuffers();
        std::swap(VBO, other.VBO);
        std::swap(VAO, other.VAO);
        std::swap(EBO, other.EBO);
        std::swap(vertices, other.vertices);
        std::swap(indices, other.indices);
    }
    return *this;
}

// 800 622 1147
void Mesh::deleteBuffers(){
    glDeleteVertexArrays(1, &VAO);
    glDeleteBuffers(1, &VBO);
    glDeleteBuffers(1, &EBO);
    VBO = VAO = EBO = 0;
    vertices.clear();
    indices.clear();
}

void Mesh::calculateVertices(){
    vertices = {
        // positions         // colors          // texcoords

         0.5f,  0.5f, 0.5f,  0.0f, 0.0f, 0.0f,  1.0, 0.0, // top right
         0.5f, -0.5f, 0.5f,  1.0f, 0.0f, 0.0f,  1.0, 1.0, // bottom right
        -0.5f, -0.5f, 0.5f,  0.0f, 1.0f, 0.0f,  0.0, 1.0, // bottom left
        -0.5f,  0.5f, 0.5f,  1.0f, 0.0f, 1.0f,  0.0, 0.0 // top left
    };

    indices = {  // note that we start from 0!
        0, 1, 3,  // first Triangle
        1, 2, 3   // second Triangle
    };

}

void Mesh::initializeBuffers(){
    // bind the Vertex Array Object first, then bind and set vertex buffer(s), and then configure vertex attributes(s).

    glGenVertexArrays(1, &VAO);
    glGenBuffers(1, &VBO);
    glGenBuffers(1, &EBO);

    glBindVertexArray(VAO);

    glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(),  GL_STATIC_DRAW);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),  GL_STATIC_DRAW);

    // position attribute
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    glEnableVertexAttribArray(0);

    // color attribute
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    glEnableVertexAttribArray(1);

    // texcoord attribute
    glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    glEnableVertexAttribArray(2);

}

void Mesh::render(){ // shouldnt be here?
    // glBindVertexArray(VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
    // glDrawArrays(GL_TRIANGLES, 0, 3);
    glDrawElements(GL_TRIANGLES, indices.size(), GL_UNSIGNED_INT, 0);
}


void Mesh::printVertices(){
    std::cout << "      Position                 Color                   TexCoord" << std::endl;
    for (int i = 0; i < vertices.size();){
        std::printf("%3d  ", i/8);
        for (int j = 0; j < 8; j++, i++){
            std::printf(" %6.3f ",  vertices[i]);
        }
        std::cout << std::endl;
    }

    std::cout << "\nIndices" << std::endl;
    for (int i = 0; i < indices.size();){
        std::printf("%3d  ", i/3);
        for (int j = 0; j < 3; j++, i++){
            std::printf(" %d ", indices[i]);
        }
        std::cout << std::endl;
    }

}

