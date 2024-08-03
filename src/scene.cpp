#include <scene.hpp>

#include <shader_s.hpp>
#include <mesh.hpp>
#include <camera.hpp>

Scene::Scene(){
    meshes = std::vector<Mesh>();
    shaders = std::vector<Shader>();
    cameras = std::vector<Camera>();
}

Scene::~Scene(){
    // DeleteSceneObjects();
}

void Scene::DeleteSceneObjects(){
    meshes.clear();
    shaders.clear();
    cameras.clear();
}

Scene::Scene(Scene&& other) : meshes(std::move(other.meshes)), shaders(std::move(other.shaders)), cameras(std::move(other.cameras)){
    other.meshes.clear();
    other.shaders.clear();
    other.cameras.clear();
}

Scene& Scene::operator=(Scene&& other){
    if (this != &other){
        DeleteSceneObjects();
        meshes = std::move(other.meshes);
        shaders = std::move(other.shaders);
        cameras = std::move(other.cameras);
        other.meshes.clear();
        other.shaders.clear();
        other.cameras.clear();
    }
    return *this;
}

// bool Scene::RegisterEntity(Mesh* mesh, Shader* shader){
//     this->meshes.push_back(mesh);
//     this->shaders.push_back(shader);

//     return true;
// }


// bool Scene::RegisterCamera(Camera * camera){
//     this->cameras.push_back(camera);
//     return true;
// }


Scene Scene::GenerateDefaultScene(){
    Scene retScene;

    retScene.meshes.push_back(Mesh::CreateCube());
    retScene.shaders.push_back(Shader("./resources/shader/base.vert", "./resources/shader/base.frag"));
    retScene.cameras.push_back(Camera(800.0, 600.0));
    // Camera camera((float) SCR_WIDTH, (float) SCR_HEIGHT);


    // // Need to make these RAII and figure out how to provide inputs
    // Shader ourShader("./resources/shader/exv.vert", "./resources/shader/exf.frag");
    // // Shader baseShader("./resources/shader/base.vert", "./resources/shader/base.frag");
    // Shader seaShader("./resources/shader/sea.vert", "./resources/shader/sea.frag");
    // Shader seafloorShader("./resources/shader/seafloor.vert", "./resources/shader/seafloor.frag");
    // // assert(false);

    // ourShader.use();
    // ourShader.setInt("texture1", 0);
    // ourShader.setInt("texture2", 1);

    // seaShader.use();
    // seaShader.setInt("texture1", 0);
    // seaShader.setInt("texture2", 1);



    // glm::mat4 view = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    // glm::mat4 projection = camera.getProjMatrix();

    // // // camera/view transformation
    // Cube skybox;
    // skybox.scale(glm::vec3(100.0));

    // Cube xcoord;
    // xcoord.translate(glm::vec3(10.0, 0.0, 0.0));

    // Cube ycoord;
    // ycoord.translate(glm::vec3(0.0, 10.0, 0.0));

    // Cube zcoord;
    // zcoord.translate(glm::vec3(0.0, 0.0, 10.0));


    // Cube light;
    // glm::vec3 lightPos = glm::vec3(0.0, 15.0, 5.0);
    // light.translate(lightPos);

    // Plane sea(glm::ivec3(1024));
    // sea.scale(glm::vec3(100.0));

    // Plane floor(glm::ivec2(1024));
    // floor.scale(glm::vec3(100.0, 50.0, 100.0));
    // floor.translate(glm::vec3(0.0, -0.5, 0.0));

    // retScene.meshes.push_back({skybox, xcoord, ycoord, zcoord, light, sea, floor}); // uses move semantics
    // retScene.cameras.push_back(camera); // literally just data
    // retScene.shaders.push_back(ourShader);


    return retScene;
}