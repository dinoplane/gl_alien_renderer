#include <macros.h>
#include <mesh.hpp>

#include <iostream>
void Mesh::GenerateBuffers(Mesh* mesh, const std::vector<Vertex> &vertices, const std::vector<uint> &indices){
    // Flexible because I can actually use this to batch the buffer creation!

    glCreateVertexArrays(1, &mesh->VAO);
    // glGenVertexArrays(1, &VAO);
    // glBindVertexArray(VAO);

    glCreateBuffers(1, &mesh->VBO);
    // glGenBuffers(1, &VBO);
    // glBindBuffer(GL_ARRAY_BUFFER, VBO);
    glNamedBufferStorage(mesh->VBO, vertices.size() * sizeof(Vertex), vertices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(),  GL_STATIC_DRAW);

    glCreateBuffers(1, &mesh->EBO);
    // glGenBuffers(1, &EBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
    glNamedBufferStorage(mesh->EBO, indices.size() * sizeof(uint), indices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),  GL_STATIC_DRAW);

    glVertexArrayVertexBuffer(mesh->VAO, 0, mesh->VBO, 0, sizeof(Vertex));
    glVertexArrayElementBuffer(mesh->VAO, mesh->EBO);

    glEnableVertexArrayAttrib(mesh->VAO, 0);
    glEnableVertexArrayAttrib(mesh->VAO, 1);
    glEnableVertexArrayAttrib(mesh->VAO, 2);

    // position attribute
    glVertexArrayAttribFormat(mesh->VAO, 0, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, position));
    glVertexArrayAttribBinding(mesh->VAO, 0, 0);
    // glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)0);
    // glEnableVertexAttribArray(0);

    std::cout << "Position offset: " << offsetof(Vertex, position) << std::endl;
    // normal attribute
    glVertexArrayAttribFormat(mesh->VAO, 1, 3, GL_FLOAT, GL_FALSE, offsetof(Vertex, normal));
    glVertexArrayAttribBinding(mesh->VAO, 1, 0);

    std::cout << "Normal offset: " << offsetof(Vertex, normal) << std::endl;
    // glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(3 * sizeof(float)));
    // glEnableVertexAttribArray(1);

    // texcoord attribute
    glVertexArrayAttribFormat(mesh->VAO, 2, 2, GL_FLOAT, GL_FALSE, offsetof(Vertex, texcoords));
    glVertexArrayAttribBinding(mesh->VAO, 2, 0);
    std::cout << "Texcoords offset: " << offsetof(Vertex, texcoords) << std::endl;
    // glVertexAttribPointer(2, 2, GL_FLOAT, GL_FALSE, 8 * sizeof(float), (void*)(6 * sizeof(float)));
    // glEnableVertexAttribArray(2);

    std::cout << "sizeof(Vertex) " << sizeof(Vertex) << std::endl;
    std::cout << "sizeof(vertices.data()) " << sizeof(vertices.data()) << std::endl;
    std::cout << "sizeof(vertices.data[0]) " << sizeof(vertices.data()[0]) << std::endl;
}

Mesh Mesh::CreateCube(){
    const std::vector<Vertex> vertices ({
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
    });
std::vector<Vertex> verticesVec (vertices);
    const std::vector<uint> indices {  // note that we start from 0!
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
    };

    Mesh mesh;
    GenerateBuffers(&mesh, verticesVec, indices);
    mesh.indexCount = indices.size();
    return mesh;
}
