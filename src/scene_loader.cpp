#include <scene_loader.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include <camera.hpp>
#include <volume.hpp>

static glm::vec3 ParseVec3(const std::string& vec3Str)
{
    glm::vec3 vec;
    std::istringstream iss(vec3Str);
    iss >> vec.x >> vec.y >> vec.z;
    return vec;
}

Scene SceneLoader::LoadScene(const SceneData& sceneData)
{
    Scene scene;
    // Add entities as instances. 
    for (const EntityData& entityData : sceneData.entitiesData)
    {
        std::string classname = entityData.className;
        if ( classname == "camera" ){
            Camera camera;
            glm::vec3 cameraAngles = Pa
            camera.fovY = std::stof(entityData.kvps.at("fovy"));
            camera.zNear = std::stof(entityData.kvps.at("near"));
            camera.zFar = std::stof(entityData.kvps.at("far"));
            camera.position = ParseVec3(entityData.kvps.at("origin"));
            rseVec3(entityData.kvps.at("angles"));
            camera.pitch = cameraAngles.x;
            camera.yaw = cameraAngles.y;
            scene.initCamConfigs.push_back(camera);

        } else {
            EntityInstanceDatax sceneEntity = scene.entityInstanceMap[];
            // if (entityData != scene.entityInstanceMap.end()){
            //     // Add instance
            sceneEntity.modelToWorldMat.push_back(entityData.transform.GetModelMatrix());
            sceneEntity.instMesh = Mesh::CreateCube();

            sceneEntity.boundingVolumes.push_back(sceneEntity.instMesh.boundingVolume->ToGPUSphere()); // This does not spark joy, the bounding volume should just be part of the mesh :skull: 
            sceneEntity.isInstMeshRendered.push_back(1);
            sceneEntity.instCount += 1;

            
            
            // }
        }
    //     Entity entity;
    //     entity.mesh = new Mesh(entityData.meshData.vertices, entityData.meshData.indices);
    //     entity.transform = entityData.transform;
    //     scene.entities.push_back(entity);
    }

    // for (const Camera& camera : sceneData.cameraData)
    // {
    //     scene.initCamConfigs.push_back(camera);
    // }

    return scene;
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
            } else {

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
        //     iss >> front.x >> front.y >> front.z;
        //     camera = Camera(fov, aspect, near, far, pos, up, front);
        //     sceneData.cameraData.push_back(camera);
        // }
        }

        }


    }
    


    PrintSceneData(sceneData);
    return sceneData;
}