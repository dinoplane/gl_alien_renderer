#include <mesh.hpp>
#include <util.h>
#include <camera.hpp>

#include <iostream>

Mesh Mesh::CreateCube(){
    Mesh mesh;
    mesh.primitives.push_back(Primitive::CreateCube());
    mesh.boundingVolume = Sphere(glm::vec3(0.f, 0.f, 0.f), 0.866f);
    return mesh;
}

Mesh Mesh::CreateFrustum(const Camera& cam){
    Mesh mesh;
    mesh.primitives.push_back(Primitive::CreateFrustum(cam));
    mesh.boundingVolume = Sphere(glm::vec3(0.f, 0.f, 0.f), 0.866f);
    
    return mesh;
}

Mesh Mesh::CreatePyramid(){
    Mesh mesh;
    mesh.primitives.push_back(Primitive::CreatePyramid());
    mesh.boundingVolume = Sphere(glm::vec3(0.f, 0.f, 0.f), 0.707f);
    return mesh;
}


void Primitive::GenerateBuffers(Primitive* primitive, const std::vector<Vertex> &vertices, const std::vector<uint32_t> &indices){
    primitive->indexCount = indices.size();
    glCreateBuffers(1, &primitive->VBO);
    // glGenBuffers(1, &primitive->VBO);
    // glBindBuffer(GL_ARRAY_BUFFER, primitive->VBO);
    glNamedBufferStorage(primitive->VBO, vertices.size() * sizeof(Vertex), vertices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(),  GL_STATIC_DRAW);

    glCreateBuffers(1, &primitive->EBO);
    // glGenBuffers(1, &primitive->EBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, primitive->EBO);
    glNamedBufferStorage(primitive->EBO, indices.size() * sizeof(uint32_t), indices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),  GL_STATIC_DRAW);
}

void Primitive::GenerateDebugBuffers(Primitive* primitive, const std::vector<glm::vec3> &vertices, const std::vector<uint32_t> &indices){

    primitive->indexCount = indices.size();
    glCreateBuffers(1, &primitive->VBO);
    // glGenBuffers(1, &primitive->VBO);
    // glBindBuffer(GL_ARRAY_BUFFER, primitive->VBO);
    glNamedBufferStorage(primitive->VBO, vertices.size() * sizeof(glm::vec3), vertices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(),  GL_STATIC_DRAW);

    glCreateBuffers(1, &primitive->EBO);
    // glGenBuffers(1, &primitive->EBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, primitive->EBO);
    glNamedBufferStorage(primitive->EBO, indices.size() * sizeof(uint32_t), indices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),  GL_STATIC_DRAW);
}

// void Mesh::RebindDebug(Mesh* primitive){
//     glBindVertexArray(mesh->VAO);
//     glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
//     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);

//     glEnableVertexArrayAttrib(mesh->VAO, 0);

//     // position attribute
//     glVertexArrayAttribFormat(mesh->VAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
//     glVertexArrayAttribBinding(mesh->VAO, 0, 0);
// }

Primitive Primitive::CreateCube(){
    const std::vector<Vertex> vertices ({ // counterclockwise
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
    const std::vector<uint32_t> indices {  // note that we start from 0!
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

    Primitive mesh;
    GenerateBuffers(&mesh, vertices, indices);
    
    return mesh;
}



Primitive Primitive::CreateFrustum(const Camera& cam){
    const float aspect = cam.width / cam.height;
    const float halfVSide = cam.zFar * tanf(glm::radians(cam.fovY) * .5f);
    const float halfHSide = halfVSide * aspect;
    const float nearToFarRatio = cam.zNear / cam.zFar;

    std::vector<glm::vec3> vertices({
        { 0.0,  0.0,  0.0},
        {cam.zFar,  halfVSide,  halfHSide},
        {cam.zFar,  halfVSide, -halfHSide},
        {cam.zFar, -halfVSide, -halfHSide},
        {cam.zFar, -halfVSide,  halfHSide},
    });

    const std::vector<uint32_t> indices({
        // Front face (near plane)
        1,  3, 4,
        3,  1, 2,  // First square
        // 1, 3, 2, 4, // Diagonals
        0,  3, 2,        // Left face
        0,  1, 4,        // Right face
        0,  2, 1,        // Top face
        0,  4, 3,         // Bottom face
        });

    Primitive mesh;
    GenerateDebugBuffers(&mesh, vertices, indices);
    return mesh;
}

Primitive Primitive::CreatePyramid(){
    const std::vector<Vertex>
    vertices ({
        // positions                // colors              // texcoords
        {{ 0.5f, -0.5f,  0.5f},    { 0.0f, 1.0f,  1.0f},   {1.0f, 1.0f}},
        {{-0.5f, -0.5f,  0.5f},    { 0.0f, 1.0f,  1.0f},   {0.0f, 1.0f}},
        {{ 0.0f,  0.5f,  0.0f},    { 0.0f, 1.0f,  1.0f},   {0.5f, 0.0f}},

        {{-0.5f, -0.5f,  0.5f},    {-1.0f, 1.0f,  0.0f},   {0.0f, 1.0f}},
        {{-0.5f, -0.5f, -0.5f},    {-1.0f, 1.0f,  0.0f},   {1.0f, 1.0f}},
        {{ 0.0f,  0.5f,  0.0f},    {-1.0f, 1.0f,  0.0f},   {0.5f, 0.0f}},

        {{-0.5f, -0.5f, -0.5f},    { 0.0f, 1.0f, -1.0f},   {1.0f, 1.0f}},
        {{ 0.5f, -0.5f, -0.5f},    { 0.0f, 1.0f, -1.0f},   {0.0f, 1.0f}},
        {{ 0.0f,  0.5f,  0.0f},    { 0.0f, 1.0f, -1.0f},   {0.5f, 0.0f}},

        {{ 0.5f, -0.5f, -0.5f},    { 1.0f, 1.0f,  0.0f},   {0.0f, 1.0f}},
        {{ 0.5f, -0.5f,  0.5f},    { 1.0f, 1.0f,  0.0f},   {1.0f, 1.0f}},
        {{ 0.0f,  0.5f,  0.0f},    { 1.0f, 1.0f,  0.0f},   {0.5f, 0.0f}},

        {{ 0.5f, -0.5f,  0.5f},    {-1.0f, 0.0f,  0.0f},   {0.0f, 1.0f}},
        {{-0.5f, -0.5f,  0.5f},    {-1.0f, 0.0f,  0.0f},   {1.0f, 1.0f}},
        {{-0.5f, -0.5f, -0.5f},    {-1.0f, 0.0f,  0.0f},   {0.0f, 1.0f}},
        {{ 0.5f, -0.5f, -0.5f},    {-1.0f, 0.0f,  0.0f},   {1.0f, 1.0f}},
    });

    const std::vector<uint32_t> indices ({  // note that we start from 0!
        0, 1, 2,  // front tr br tl

        3, 4, 5,  // right tr br tl

        6, 7, 8,  // back tr br tl

        9, 10, 11,  // left tr br tl

        12, 13, 14,  // bot tr br tl
        14, 15, 12   // bot br bl tl
    });
    Primitive mesh;
    GenerateBuffers(&mesh, vertices, indices);

    return mesh;
}
