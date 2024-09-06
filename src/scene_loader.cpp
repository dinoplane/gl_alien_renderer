#include <scene_loader.hpp>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <glm/gtc/type_ptr.hpp>


#include <glm/gtx/string_cast.hpp>
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

static glm::vec3 ParseVec3(const std::string& vec3Str)
{
    glm::vec3 vec;
    std::istringstream iss(vec3Str);
    iss >> vec.x >> vec.y >> vec.z;
    return vec;
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
            data.kvps.push_back({key, value});
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