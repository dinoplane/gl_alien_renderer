#ifndef SCENE_DATA_H
#define SCENE_DATA_H

#include <util.h>
#include <vector>
#include <string>
#include <transform.hpp>


class Camera;

struct MeshData {
    std::vector<Vertex> vertices;
    std::vector<unsigned int> indices;

    static MeshData CreateCube();
};

struct EntityData {
    MeshData meshData;
    Transform transform;
};

struct SceneData {
    std::vector<EntityData> entitiesData;
    std::vector<std::pair<std::string, std::string>> shaderPaths;
    std::vector<Camera> cameraData;


    static SceneData GenerateDefaultScene();

    // at this point its becoming an engine LMAO
    // static SceneData GenerateSceneFromConfig(const std::string& configPath);
    // If I process the file multithreadedly, TECHNICALLY, i can guarantee a random bvh insertion...
    // Wait thats lk troll and funny.
};

#endif
