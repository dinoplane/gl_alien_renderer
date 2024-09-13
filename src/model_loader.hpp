#ifndef MODEL_LOADER_H
#define MODEL_LOADER_H
#include <filesystem>


#include <mesh.hpp>

#include <fastgltf/core.hpp>
#include <fastgltf/types.hpp>
#include <fastgltf/tools.hpp>

class ModelLoader {
    public: 
    static bool LoadGLTF(const std::filesystem::path path, fastgltf::Asset* retAsset);

    static bool LoadMesh(const fastgltf::Asset& asset, const fastgltf::Mesh& mesh, Mesh* outMesh);
};

#endif