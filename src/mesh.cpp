#include <util.h>
#include <mesh.hpp>
#include <camera.hpp>

#include <iostream>


void Mesh::GenerateBuffers(Mesh* mesh, const std::vector<Vertex> &vertices, const std::vector<uint> &indices){
    glCreateBuffers(1, &mesh->VBO);
    // glGenBuffers(1, &mesh->VBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    glNamedBufferStorage(mesh->VBO, vertices.size() * sizeof(Vertex), vertices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(Vertex), vertices.data(),  GL_STATIC_DRAW);

    glCreateBuffers(1, &mesh->EBO);
    // glGenBuffers(1, &mesh->EBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
    glNamedBufferStorage(mesh->EBO, indices.size() * sizeof(uint), indices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),  GL_STATIC_DRAW);
}

void Mesh::GenerateDebugBuffers(Mesh* mesh, const std::vector<glm::vec3> &vertices, const std::vector<uint> &indices){
    glCreateBuffers(1, &mesh->VBO);
    // glGenBuffers(1, &mesh->VBO);
    // glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
    glNamedBufferStorage(mesh->VBO, vertices.size() * sizeof(glm::vec3), vertices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ARRAY_BUFFER, vertices.size() * sizeof(float), vertices.data(),  GL_STATIC_DRAW);

    glCreateBuffers(1, &mesh->EBO);
    // glGenBuffers(1, &mesh->EBO);
    // glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);
    glNamedBufferStorage(mesh->EBO, indices.size() * sizeof(uint), indices.data(), GL_DYNAMIC_STORAGE_BIT);
    // glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), indices.data(),  GL_STATIC_DRAW);
}

// void Mesh::RebindDebug(Mesh* mesh){
//     glBindVertexArray(mesh->VAO);
//     glBindBuffer(GL_ARRAY_BUFFER, mesh->VBO);
//     glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, mesh->EBO);

//     glEnableVertexArrayAttrib(mesh->VAO, 0);

//     // position attribute
//     glVertexArrayAttribFormat(mesh->VAO, 0, 3, GL_FLOAT, GL_FALSE, 0);
//     glVertexArrayAttribBinding(mesh->VAO, 0, 0);
// }

Mesh Mesh::CreateCube(){
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
    GenerateBuffers(&mesh, vertices, indices);
    mesh.indexCount = indices.size();
    return mesh;
}



Mesh Mesh::CreateFrustum(const Camera& cam){


    const float aspect = cam.width / cam.height;
    const float halfVSide = cam.zFar * tanf(glm::radians(cam.fovY) * .5f);
    const float halfHSide = halfVSide * aspect;
    const float nearToFarRatio = cam.zNear / cam.zFar;
    // std::vector<glm::vec3> vertices({
    //     { halfHSide,  halfVSide, cam.zNear},
    //     {-halfHSide,  halfVSide, cam.zNear},
    //     {-halfHSide, -halfVSide, cam.zNear},
    //     { halfHSide, -halfVSide, cam.zNear},

    //     { halfHSide,  halfVSide, cam.zFar*0.1},
    //     {-halfHSide,  halfVSide, cam.zFar*0.1},
    //     {-halfHSide, -halfVSide, cam.zFar*0.1},
    //     { halfHSide, -halfVSide, cam.zFar*0.1},
    // });



    // std::vector<glm::vec3> vertices({
    //     { 0.5,  0.5, -0.5},
    //     {-0.5,  0.5, -0.5},
    //     {-0.5, -0.5, -0.5},
    //     { 0.5, -0.5, -0.5},

    //     { 0.5,  0.5, 0.5},
    //     {-0.5,  0.5, 0.5},
    //     {-0.5, -0.5, 0.5},
    //     { 0.5, -0.5, 0.5},
    // });


    std::vector<glm::vec3> vertices({
        { 0.0,  0.0,  0.0},
        {cam.zFar,  halfVSide,  halfHSide},
        {cam.zFar,  halfVSide, -halfHSide},
        {cam.zFar, -halfVSide, -halfHSide},
        {cam.zFar, -halfVSide,  halfHSide},
    });

    // std::vector<glm::vec3> vertices({
    //     { 0.0,  0.0,  0.0},
    //     {5.0,  1.5,   1.5},
    //     {5.0,  1.5,  -1.5},
    //     {5.0, -1.5,  -1.5},
    //     {5.0, -1.5,   1.5},
    // });


    std::cout << halfHSide << " " << halfVSide << " " << cam.zNear << " " << cam.zFar << std::endl;
    /*
    // for (int i = 0; i < 4; ++i){
    //     vertices.push_back({halfHSide, halfVSide, cam.zNear});

    //     currMult = currMult >> 1;
    // }

    // uint currMult = 12;
    // for (int i = 0; i < 4; ++i){
    //     vertices.push_back({halfHSide, halfVSide, cam.zNear});

    //     currMult = currMult >> 1;
    // }
    */
    /*
        std::vector<uint> indices({
            // Front face (near plane)
            0, 1, 2,  // First triangle
            0, 2, 3,  // Second triangle

            // Back face (far plane)
            5, 4, 7,  // First triangle
            5, 7, 6,  // Second triangle

            // Left face
            1, 5, 6,  // First triangle
            1, 6, 2,  // Second triangle

            // Right face
            4, 0, 3,  // First triangle
            4, 3, 7,  // Second triangle

            // Top face
            4, 5, 1,  // First triangle
            4, 1, 0,  // Second triangle

            // Bottom face
            3, 6, 2,  // First triangle
            3, 7, 6   // Second triangle
        });

        std::vector<uint> indices({
            // Front face (near plane)
            0, 1, 1, 2, 2, 3, 3, 0,  // First square
            1, 5, 5, 6, 6, 2, 2, 1,  // Second square
            5, 4, 4, 7, 7, 6, 6, 5,  // Third square
            4, 0, 0, 3, 3, 7, 7, 4,  // Fourth square
            4, 5, 5, 1, 1, 0, 0, 4,  // Fifth square
            3, 2, 2, 6, 6, 7, 7, 3   // Sixth square
        });
    */


    // std::vector<uint> indices({
    //     // Front face (near plane)
    //     1, 4, 4, 3, 3, 2, 2, 1,  // First square
    //     1, 3, 2, 4, // Diagonals
    //     0, 2, 2, 3, 3, 0,        // Left face
    //     0, 4, 4, 1, 1, 0,        // Right face
    //     0, 1, 1, 2, 2, 0,        // Top face
    //     0, 3, 3, 4, 4, 0         // Bottom face
    //     });

    const std::vector<uint> indices({
        // Front face (near plane)
        1,  3, 4,
        3,  1, 2,  // First square
        // 1, 3, 2, 4, // Diagonals
        0,  3, 2,        // Left face
        0,  1, 4,        // Right face
        0,  2, 1,        // Top face
        0,  4, 3,         // Bottom face
        });



    Mesh mesh;
    GenerateDebugBuffers(&mesh, vertices, indices);
    mesh.indexCount = indices.size();
    mesh.boundingVolume = new Sphere(glm::vec3(0.f, 0.f, 0.f), 0.707f);
    // mesh.boundingVolume =
    return mesh;
}

Mesh Mesh::CreatePyramid(){
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



    const std::vector<uint> indices ({  // note that we start from 0!
        0, 1, 2,  // front tr br tl

        3, 4, 5,  // right tr br tl

        6, 7, 8,  // back tr br tl

        9, 10, 11,  // left tr br tl

        12, 13, 14,  // bot tr br tl
        14, 15, 12   // bot br bl tl
    });
    Mesh mesh;
    GenerateBuffers(&mesh, vertices, indices);
    mesh.indexCount = indices.size();

    mesh.boundingVolume = new Sphere(glm::vec3(0.f, 0.f, 0.f), 0.86f);
    return mesh;
}
