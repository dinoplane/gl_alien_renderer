#ifndef MODEL_LOADER_H
#define MODEL_LOADER_H
#include <filesystem>

#include <fastgltf/core.hpp>
#include <fastgltf/types.hpp>
#include <fastgltf/tools.hpp>



struct Mesh;
class Model;
struct Material;
struct Texture;

class ModelLoader {

    public: 

        static bool LoadGLTF(const std::filesystem::path path, fastgltf::Asset* retAsset);

        static bool LoadMesh(const fastgltf::Asset& asset, const fastgltf::Mesh& mesh, const Model& model, Mesh* outMesh);

        static bool LoadModel(const fastgltf::Asset& asset, Model *model);

        // static bool LoadTexture(const fastgltf::Asset& asset, const fastgltf::Texture& texture);

        static bool LoadImage(const fastgltf::Asset& asset, const fastgltf::Image& image, Texture* outTexture);

        static bool LoadMaterial(const fastgltf::Asset& asset, const fastgltf::Material& material, Material* outMaterial);

};

#endif