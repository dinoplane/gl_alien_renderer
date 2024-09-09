#ifndef SCENE_DATA_H
#define SCENE_DATA_H

#include <util.h>
#include <vector>
#include <string>
#include <unordered_map>

#include <transform.hpp>

#include <camera.hpp>

struct EntityData {
    
    std::string className;
    //MeshData meshData;
    Transform transform;
    bool isInstance;
    std::unordered_map<std::string, std::string> kvps;
};

struct SceneData {
    // std::set<std::string> entityKeys; // so i could conserve space here but uh... save for another day 
    std::vector<EntityData> entitiesData;
    std::vector<std::pair<std::string, std::string>> shaderPaths;
    std::vector<Camera> cameraData;

    bool errFlag = false;

    // at this point its becoming an engine LMAO
    // static SceneData GenerateSceneFromConfig(const std::string& configPath);
    // If I process the file multithreadedly, TECHNICALLY, i can guarantee a random bvh insertion...
    // Wait thats lk troll and funny.
};


void PrintSceneData(const SceneData& sceneData);

#endif
