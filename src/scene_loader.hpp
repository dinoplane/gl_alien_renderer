#ifndef SCENE_LOADER_H
#define SCENE_LOADER_H

#include <vector>
#include <string>
#include <unordered_map>
#include <scene.hpp>
#include <scenedata.hpp>


class SceneLoader {
    public: 
    static Scene LoadScene(const SceneData& sceneData);
    static SceneData LoadSceneData(const std::string& scenePath);


};



#endif