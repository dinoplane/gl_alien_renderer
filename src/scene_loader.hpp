#ifndef SCENE_LOADER_H
#define SCENE_LOADER_H

#include <vector>
#include <string>
#include <unordered_map>
#include <scene.hpp>
#include <scenedata.hpp>


class SceneLoader {
    public: 
    static void LoadScene(const SceneData& sceneData, Scene* scene);
    static SceneData LoadSceneData(const std::string& scenePath);
    static void LoadCameraSettings(const std::string& cameraSettingsPath, std::vector<Camera>* camerasm, std::vector<Primitive>* debugMeshes);
    static void SaveCameraSettings(const std::string& cameraSettingsPath, const std::vector<Camera>& cameras);
};



#endif