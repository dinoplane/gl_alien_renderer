#include <scene.hpp>

#include <shader_s.hpp>
#include <mesh.hpp>
#include <transform.hpp>
#include <entity.hpp>
#include <camera.hpp>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <iostream>
#include <glm/gtx/string_cast.hpp>
#include <gpu_structs.hpp>


void EntityInstanceData::GenerateInstanceBuffers(){ // TODO honestly i could make this function take in the vectors.
    glCreateBuffers(1, &instModelMatrixBuffer);
    glNamedBufferStorage(instModelMatrixBuffer, modelToWorldMat.size() * sizeof(glm::mat4), (const void *) modelToWorldMat.data(), GL_DYNAMIC_STORAGE_BIT);

    // glCreateBuffers(1, &instBoundingVolumeBuffer);
    // glNamedBufferStorage(instBoundingVolumeBuffer, boundingVolumes.size() * sizeof(GPUSphere), (const void *) boundingVolumes.data(), GL_DYNAMIC_STORAGE_BIT);

    isInstMeshRendered.resize(instCount, 1); // TODO theres gotta be a way to initialize this on the GPU
    // fmt::print("Size of isInstMeshRendered: {}\n", sizeof(GLboolean) * );
    fmt::print("Size of Glboolean: {}\n", sizeof(GLboolean));

    glCreateBuffers(1, &instMeshRenderedBuffer);
    glNamedBufferStorage(instMeshRenderedBuffer, isInstMeshRendered.size() * sizeof(int), (const void *) isInstMeshRendered.data(), GL_DYNAMIC_STORAGE_BIT);

    GLint maxComputeWorkGroupCount;
    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &maxComputeWorkGroupCount);
    fmt::print("Max Compute Work Group Count x: {}\n", maxComputeWorkGroupCount);
}

Scene::Scene(){
    entities = std::vector<Entity>();
    shaders = std::vector<Shader>();
    initCamConfigs = std::vector<Camera>();
    entityInstanceMap = std::unordered_map<std::string, EntityInstanceData>();
}

Scene::~Scene(){
    // DeleteSceneObjects();
}

void Scene::DeleteSceneObjects(){
    entities.clear();
    shaders.clear();
    initCamConfigs.clear();
    debugMeshes.clear();
    entityInstanceMap.clear();
}

Scene::Scene(Scene&& other) :
            entities(std::move(other.entities)),
            shaders(std::move(other.shaders)),
            initCamConfigs(std::move(other.initCamConfigs)),
            entityInstanceMap(std::move(other.entityInstanceMap))
            {
    other.entities.clear();
    other.shaders.clear();
    other.initCamConfigs.clear();
    other.debugMeshes.clear();
    other.entityInstanceMap.clear();
}

Scene& Scene::operator=(Scene&& other){
    if (this != &other){
        DeleteSceneObjects();
        entities = std::move(other.entities);
        shaders = std::move(other.shaders);
        initCamConfigs = std::move(other.initCamConfigs);
        debugMeshes = std::move(other.debugMeshes);
        entityInstanceMap = std::move(other.entityInstanceMap);
        other.entities.clear();
        other.shaders.clear();
        other.initCamConfigs.clear();
        other.debugMeshes.clear();
        other.entityInstanceMap.clear();
    }
    return *this;
}

// bool Scene::RegisterEntity(Entity* entity, Shader* shader){
//     this->entities.push_back(entity);
//     this->shaders.push_back(shader);

//     return true;
// }


// bool Scene::RegisterCamera(Camera * camera){
//     this->initCamConfigs.push_back(camera);
//     return true;
// }

void Scene::RebindAllMeshes(){
    // for (Entity& entity : entities){
    //     // Mesh::Rebind(&entity.mesh);
    // }

    // for (Mesh& mesh : debugMeshes){
    //     Mesh::RebindDebug(&mesh);
    // }
}

Scene Scene::GenerateBasicScene(){
    Scene retScene;
    Transform transform;
    transform.SetPosition(glm::vec3(0.0, 0.0, 0.0));
    transform.SetRotation(glm::vec3(0.0, 0.0, 0.0));
    transform.SetScale(glm::vec3(1.0, 1.0, 1.0));

    EntityInstanceData entityInstanceData;
    entityInstanceData.instModel = Model::CreateCube();
    entityInstanceData.instCount = 0;
    retScene.entityInstanceMap["hi"] =  entityInstanceData;
    retScene.entityInstanceMap["hi"].modelToWorldMat.push_back(transform.GetModelMatrix());
    retScene.entityInstanceMap["hi"].instCount += 1;
    retScene.entityInstanceMap["hi"].GenerateInstanceBuffers();


    retScene.shaders.push_back(Shader("./resources/shader/base.vert", "./resources/shader/base.frag"));

    // retScene.initCamConfigs.push_back(Camera(400.0, 400.0, glm::vec3(0.0, 2.0, 0.0), glm::vec3(0.0, 1.0, 0.0), -320.0f, -70.0f));
    // retScene.initCamConfigs.push_back(Camera(400.0, 400.0, glm::vec3(0.0, 2.0, 0.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f));

    return retScene;
}

// void AddEntity(Entity entity){
//     retScene.entities
// }

Scene Scene::GenerateDefaultScene(){
    const int LATTICE_LENGTH = 10;
    Scene retScene;

    Transform transform;
    EntityInstanceData entityInstanceData;
    entityInstanceData.instModel = Model::CreateCube();
    entityInstanceData.instCount = 0;
    std::cout << (entityInstanceData.instModel.boundingVolume.center.x) << std::endl;
    retScene.entityInstanceMap["hi"] =  entityInstanceData;

    transform.SetRotation(glm::vec3(0.0, 0.0, 0.0));
    transform.SetScale(glm::vec3(1.0, 1.0, 1.0));
    for (int i = 0; i < LATTICE_LENGTH; ++i){
        for (int j = 0; j < LATTICE_LENGTH; ++j){
            for (int k = 0; k < LATTICE_LENGTH; ++k){
                transform.SetPosition(glm::vec3(i * 2.0, j * 2.0, k * 2.0));
                // retScene.entities.push_back({Mesh::CreateCube(), transform});
                retScene.entityInstanceMap["hi"].modelToWorldMat.push_back(transform.GetModelMatrix());
                // retScene.entityInstanceMap["hi"].boundingVolumes.push_back(entityInstanceData.instModel.boundingVolume->ToGPUSphere());
                retScene.entityInstanceMap["hi"].instCount += 1;
            }
        }
    }
    retScene.entityInstanceMap["hi"].GenerateInstanceBuffers();


    // transform.SetRotation(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetScale(glm::vec3(1.0, 1.0, 1.0));
    // for (int i = 0; i < LATTICE_LENGTH; ++i){
    //     for (int j = 0; j < LATTICE_LENGTH; ++j){
    //         for (int k = 0; k < LATTICE_LENGTH; ++k){
    //             transform.SetPosition(glm::vec3(i * -2.0, j * -2.0, k * -2.0));
    //             retScene.entities.push_back({Mesh::CreatePyramid(), transform});
    //         }
    //     }
    // }

    // transform.SetPosition(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetRotation(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetScale(glm::vec3(1.0, 1.0, 10.0));
    // // transform.GetModelMatrix();
    // retScene.entities.push_back({Mesh::CreateCube(), transform});


    // transform.SetPosition(glm::vec3(0.0, 0.0, 0.0));
    // transform.SetRotation(glm::vec3(0.0, 90.0, 0.0));
    // transform.SetScale(glm::vec3(1.0, 10.0, 1.0));
    // retScene.entities.push_back({Mesh::CreateCube(), transform});

    retScene.shaders.push_back(Shader("./resources/shader/base.vert", "./resources/shader/base.frag"));
    retScene.shaders.push_back(Shader("./resources/shader/base_inst.vert", "./resources/shader/base_inst.frag"));
    // retScene.shaders.push_back(Shader("./resources/shader/debug.vert", "./resources/shader/debug.frag"));



    // retScene.initCamConfigs.push_back(Camera(800.0, 600.0, glm::vec3(0.0, 20.0, 0.0), glm::vec3(0.0, 1.0, 0.0), -315.0f, -60.0f));
    // retScene.initCamConfigs.push_back(Camera(800.0, 600.0));

    // for ( Camera& camera : retScene.initCamConfigs) {
    //     retScene.debugMeshes.push_back(Mesh::CreateFrustum(camera));
    // }



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

    // retScene.entities.push_back({skybox, xcoord, ycoord, zcoord, light, sea, floor}); // uses move semantics
    // retScene.initCamConfigs.push_back(camera); // literally just data
    // retScene.shaders.push_back(ourShader);


    return retScene;
}