#include <scene_loader.hpp>
#include <particle_system.hpp>
#include <cloth_system.hpp>


#include <fstream>
#include <sstream>
#include <iomanip>
#include <filesystem>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <entity.hpp>
#include <shader_s.hpp>
#include <camera.hpp>
#include <volume.hpp>
#include <model_loader.hpp>

static glm::vec3 ParseVec3(const std::string& vec3Str)
{
    glm::vec3 vec;
    std::istringstream iss(vec3Str);
    iss >> vec.x >> vec.y >> vec.z;
    return vec;
}

void SceneLoader::LoadScene(const SceneData& sceneData, Scene* scene)
{
    // Add entities as instances. 
    for (const EntityData& entityData : sceneData.entitiesData)
    {
        std::string classname = entityData.className;
        if ( classname == "camera" ){
            Camera camera;
            glm::vec3 cameraAngles = ParseVec3(entityData.kvps.at("angles"));
            camera.fovY = std::stof(entityData.kvps.at("fovy"));
            camera.zNear = std::stof(entityData.kvps.at("near"));
            camera.zFar = std::stof(entityData.kvps.at("far"));
            camera.position = ParseVec3(entityData.kvps.at("origin"));
            camera.pitch = cameraAngles.x;
            camera.yaw = cameraAngles.y;
            scene->initCamConfigs.push_back(camera);

        } else {
            const std::string meshMatKey = entityData.kvps.at("mesh") + entityData.kvps.at("material");
            if (std::stoi(entityData.kvps.at("is_instanced"))){
                auto sceneEntity = scene->entityInstanceMap.find(meshMatKey);

                if (sceneEntity == scene->entityInstanceMap.end())
                {
                    sceneEntity = scene->entityInstanceMap.insert({meshMatKey, EntityInstanceData()}).first;
                    
                    const auto& meshData = entityData.kvps.find("mesh");
                    if (meshData != entityData.kvps.end()){
                        // Load mesh
                        if (meshData->second == "cube"){
                            sceneEntity->second.instModel = Model::CreateCube();
                        } else {
                            fastgltf::Asset asset;
                            ModelLoader::LoadGLTF(std::filesystem::path(meshData->second), &asset);
                            ModelLoader::LoadModel(asset, &sceneEntity->second.instModel);
                        }
                    } else sceneEntity->second.instModel = Model::CreateCube();

                    fmt::print("Alpha cutoff: {}\n", sceneEntity->second.instModel.materials[0].alphaCutoff);
                }
                // if (entityData != scene->entityInstanceMap.end()){
                //     // Add instance
                sceneEntity->second.modelToWorldMat.push_back(entityData.transform.GetModelMatrix());
                
                // sceneEntity->second.boundingVolumes.push_back(sceneEntity->second.instMesh.boundingVolume->ToGPUSphere()); // This does not spark joy, the bounding volume should just be part of the mesh :skull: 
                sceneEntity->second.isInstMeshRendered.push_back(1);
                sceneEntity->second.instCount += 1;
            } else {
                Model mesh;
                const auto& meshData = entityData.kvps.find("mesh");
                if (meshData != entityData.kvps.end()){
                    // Load mesh
                    if (meshData->second == "cube"){
                        mesh = Model::CreateCube();
                    } else {
                        fastgltf::Asset asset;

                            ModelLoader::LoadGLTF(std::filesystem::path(meshData->second), &asset);
                            ModelLoader::LoadModel(asset, &mesh);
                    }
                } else mesh = Model::CreateCube();

                scene->entities.push_back(Entity(mesh, entityData.transform));
            }

            
            // }
        }
    //     Entity entity;
    //     entity.mesh = new Mesh(entityData.meshData.vertices, entityData.meshData.indices);
    //     entity.transform = entityData.transform;
    //     scene->entities.push_back(entity);
    }

    for (auto& [meshMatKey, entityInstanceData] : scene->entityInstanceMap){
        entityInstanceData.GenerateInstanceBuffers();
    }

    // for (const Camera& camera : sceneData.cameraData)
    // {
    //     scene->initCamConfigs.push_back(camera);
    // }
    // I could copy it over and make some edits 
    scene->shaders.push_back(Shader("./resources/shader/base.vert", "./resources/shader/base_inst.frag"));
    scene->shaders.push_back(Shader("./resources/shader/indirect_inst.vert", "./resources/shader/indirect_inst.frag"));

    fmt::printf("Size: particle systems : %u\n", scene->particleSystems.size());

    //scene->particleSystems.push_back(ParticleSystem({100, 0.0001, "gen_particle"}));
    // 
    // scene->particleSystems.push_back(ParticleSystem({100, 0.0001, "gen_particle"}));
    // ParticleSystemParameters p{ 10, 0.0001, "gen" };
    //scene->particleSystems.push_back( std::make_unique<ParticleSystem>());
    //scene->particleSystems[0]->Initialize(&p);

    // ClothSystemParameters c;
    // c.particleCount = 100;
    // c.timeStep = 0.005;
    // c.shaderName = "cloth";
    // c.clothSideLength = 20;
    // c.cellSideLength = 0.1;
    // c.totalMass = 0.01;
    // c.gravityAccel = 9.81;
    // c.youngModulus = 1.0e5;
    // c.thickness = 0.001;


    // scene->particleSystems.push_back(std::make_unique<ClothSystem>());
    // scene->particleSystems[0]->Initialize(&c);


    // scene->particleSystems.push_back();
    // ClothSystem* cs = new ClothSystem(&c);
    //ClothSystem cs2 (&c);

    // BaseParticleSystem<ClothDataBlock, ClothSystemDataBlock, ClothSystemParameters> cs2(&c);
    //BaseParticleSystem<ParticleDataBlock, ParticleSystemDataBlock, ParticleSystemParameters> cs2(&p);
    

    for (auto& ps : scene->particleSystems){
       ps->InitializeBuffers();
    }

    // scene->clothSystems.push_back(ClothSystem({100, 0.0001, "gen_particle"}));

    // for (ClothSystem& cs : scene->clothSystems){
    //    cs.InitializeBuffers();
    // }
    

    // return std::move(scene);
}


SceneData SceneLoader::LoadSceneData(const std::string& scenePath)
{
    // HONESTLY RETURN BY REFERENCE :SKULL:
    SceneData sceneData;
    // Load scene data from file
    // Here I could do a couple of things:
    // mmap the file and read it in chunks  (oo this is kinda fun to do later)
    // use istream 

    // For now, I'll just use ifstream
    std::ifstream sceneFile(scenePath);
    if (!sceneFile.is_open())
    {
        sceneData.errFlag = true;
        std::cerr << "Failed to open scene file: " << scenePath << std::endl;
        return sceneData;
    }

    std::string line;
    EntityData data;
    Camera cameraData;
    while (std::getline(sceneFile, line))
    {
        if (line.empty())
        {
            continue;
        } else if (line[0] == '{')
        {
            // Start of a new entity
            data = EntityData();
            continue;

        } else if (line[0] == '}')
        {
            // End of entity
            sceneData.entitiesData.push_back(data);
            continue;
        } else {
        std::istringstream iss(line);
        std::string key;
        std::string value;
        while (iss >> std::quoted(key))
        {
            // TODO These should be done by type not key name
            iss >> std::quoted(value);
            if (key == "classname")
            {
                data.className = value;
            } else if (key == "origin")
            {
                data.transform.SetPosition(ParseVec3(value));
            } else if (key == "angles")
            {
                data.transform.SetRotation(ParseVec3(value));
            } else if (key == "scale"){
                data.transform.SetScale(ParseVec3(value));
            }
            data.kvps[key] = value;
        // if (token == "entity")
        // {
        //     EntityData entityData;
        //     std::string meshPath;
        //     iss >> meshPath;
        //     entityData.meshData = MeshData::CreateCube();
        //     entityData.transform = Transform();
        //     sceneData.entitiesData.push_back(entityData);
        // }
        // else if (token == "camera")
        // {
        //     Camera camera;
        //     float fov, aspect, near, far;
        //     glm::vec3 pos, up, front;
        //     iss >> fov >> aspect >> near >> far;
        //     iss >> pos.x >> pos.y >> pos.z;
        //     iss >> up.x >> up.y >> up.z;
        //     iss >> front.x >> front.y >> front .z;
        //     camera = Camera(fov, aspect, near, far, pos, up, front);
        //     sceneData.cameraData.push_back(camera);
        // }
        }

        }


    }
    


    // PrintSceneData(sceneData);
    return sceneData;
}


void SceneLoader::LoadCameraSettings(const std::string& cameraSettingsPath, std::vector<Camera>* cameras, std::vector<Primitive>* debugMeshes)
{
    std::ifstream cameraSettingsFile(cameraSettingsPath);
    if (!cameraSettingsFile.is_open())
    {
        std::cerr << "Failed to open camera settings file: " << cameraSettingsPath << std::endl;
        return;
    }

    std::string line;
    while (std::getline(cameraSettingsFile, line))
    {
        if (line.empty())
        {
            continue;
        }

        std::istringstream iss(line);
        std::string token;
        //Camera camera;
        float width, height, fov, yaw, pitch, camNear, camFar;
        glm::vec3 pos, up, front;
        iss >> token >> token;
        std::cout << token << std::endl;
        iss
            >> token >> width;
        std::cout << token << std::endl;
        iss
            >> token >> height;
        iss
            >> token >> yaw;
        iss
            >> token >> pitch;
        iss
            >> token >> fov;
        iss
            >> token >> camNear;
        iss
            >> token >> camFar;
        iss >> token >> pos.x >> pos.y >> pos.z;

         cameras->push_back(Camera(width, height, pos, yaw, pitch, fov, camNear, camFar));
         debugMeshes->push_back(Primitive::CreateFrustum((*cameras)[cameras->size() - 1]));
    }

    cameraSettingsFile.close();
}

void SceneLoader::SaveCameraSettings(const std::string& cameraSettingsPath, const std::vector<Camera>& cameras)
{
    std::ofstream cameraSettingsFile(cameraSettingsPath, std::ios::out | std::ios::trunc);
    if (!cameraSettingsFile.is_open())
    {
        std::cerr << "Failed to open camera settings file: " << cameraSettingsPath << std::endl;
        return;
    }

    int cameraIdx = 0;
    for (const Camera& camera : cameras)
    {
        cameraSettingsFile  << "Camera " << cameraIdx++;
        cameraSettingsFile  << " Width: " << camera.width 
                            << " Height: " <<  camera.height 
                            << " Yaw: " << camera.yaw 
                            << " Pitch: " << camera.pitch
                            << " Fov: " << camera.fovY 
                            << " Near: " << camera.zNear 
                            << " Far: " << camera.zFar;
        cameraSettingsFile  << " Position: " << camera.position.x << " " << camera.position.y << " " << camera.position.z << std::endl;
        // cameraSettingsFile << " " << camera.front.x << " " << camera.front.y << " " << camera.front.z;
    }

    cameraSettingsFile.close();
}
