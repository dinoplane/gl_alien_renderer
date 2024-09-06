#include <scene_loader.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>

Scene SceneLoader::LoadScene(const SceneData& sceneData)
{
    Scene scene;
    // for (const EntityData& entityData : sceneData.entitiesData)
    // {
    //     Entity entity;
    //     entity.mesh = new Mesh(entityData.meshData.vertices, entityData.meshData.indices);
    //     entity.transform = entityData.transform;
    //     scene.entities.push_back(entity);
    // }

    // for (const Camera& camera : sceneData.cameraData)
    // {
    //     scene.initCamConfigs.push_back(camera);
    // }

    return scene;
}

SceneData SceneLoader::LoadSceneData(const std::string& scenePath)
{
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
    void* data;
    while (std::getline(sceneFile, line))
    {
        std::istringstream iss(line);
        std::string key;
        std::string value;
        while (iss >> std::quoted(key))
        {
            iss >> std::quoted(value)
            if (key == "classname")
            {
                EntityData entityData;
            }
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
    



    return sceneData;
}