 #include <model_loader.hpp>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"


#include <util.h>
#include <mesh.hpp>
#include <model.hpp>

#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
 
 #include <limits>
 #include <iostream>
 #include <variant>
 



bool ModelLoader::LoadGLTF(const std::filesystem::path path, fastgltf::Asset* retAsset) {
	if (!std::filesystem::exists(path)) {
		std::cout << "Failed to find " << path << '\n';
		return false;
	}

	if constexpr (std::is_same_v<std::filesystem::path::value_type, wchar_t>) {
		std::wcout << "Loading " << path << '\n';
	} else {
		std::cout << "Loading " << path << '\n';
	}

    // Parse the glTF file and get the constructed asset
    {
		static constexpr auto supportedExtensions =
			fastgltf::Extensions::KHR_mesh_quantization |
			fastgltf::Extensions::KHR_texture_transform |
			fastgltf::Extensions::KHR_materials_variants |
			fastgltf::Extensions::KHR_materials_specular;

        fastgltf::Parser parser(supportedExtensions);

        constexpr auto gltfOptions =
            fastgltf::Options::DontRequireValidAssetMember |
            fastgltf::Options::AllowDouble |
            // fastgltf::Options::LoadGLBBuffers |
            fastgltf::Options::LoadExternalBuffers |
            fastgltf::Options::LoadExternalImages |
			fastgltf::Options::GenerateMeshIndices;

		auto gltfFile = fastgltf::MappedGltfFile::FromPath(path);
		if (!bool(gltfFile)) {
			std::cout << "Failed to open glTF file: " << fastgltf::getErrorMessage(gltfFile.error()) << '\n';
			return false;
		}

        auto asset = parser.loadGltf(gltfFile.get(), path.parent_path(), gltfOptions);
        if (asset.error() != fastgltf::Error::None) {
            std::cout << "Failed to load glTF: " << fastgltf::getErrorMessage(asset.error()) << '\n';
            return false;
        }
		*retAsset = std::move(asset.get());

        }

    return true;
}

bool ModelLoader::LoadMesh(const fastgltf::Asset& asset, const fastgltf::Mesh& mesh, Model* model, Mesh* outMesh, std::vector<Vertex>* vertices, std::vector<uint>* indices) {
    // Mesh outMesh = {};
	// Mesh retMesh;
	glm::vec3 min = glm::vec3(std::numeric_limits<float>::max());
	glm::vec3 max = glm::vec3(std::numeric_limits<float>::min());

	outMesh->startPrimIdx = static_cast<uint>(model->primitives.size());
	outMesh->primCount = static_cast<uint>(mesh.primitives.size());


    for (auto it = mesh.primitives.begin(); it != mesh.primitives.end(); ++it) {
		Primitive outPrimitive;
	

		auto* positionIt = it->findAttribute("POSITION");
		assert(positionIt != it->attributes.end()); // A mesh primitive is required to hold the POSITION attribute.
		assert(it->indicesAccessor.has_value()); // We specify GenerateMeshIndices, so we should always have indices


		std::size_t baseColorTexcoordIndex = 0;


        // Get the output primitive
        auto index = std::distance(mesh.primitives.begin(), it);

		model->primitivePropertiesBufferVec.push_back( { 0, 0 } ); // I could also pack both into a value
        if (it->materialIndex.has_value()) {
            outPrimitive.materialUniformsIndex = it->materialIndex.value() + 1; // Adjust for default material
			model->primitivePropertiesBufferVec.back().materialIdx = it->materialIndex.value() + 1;

            auto& material = asset.materials[it->materialIndex.value()];

			auto& baseColorTexture = material.pbrData.baseColorTexture;
            if (baseColorTexture.has_value()) {
                auto& texture = asset.textures[baseColorTexture->textureIndex];
				if (!texture.imageIndex.has_value())
					return false;
                outPrimitive.albedoTexture = model->textures[texture.imageIndex.value()].texture;

				model->primitivePropertiesBufferVec.back().textureIdx = texture.imageIndex.value();
				// fmt::print("Albedo Loaded {}\n", outPrimitive.albedoTexture);
				if (baseColorTexture->transform && baseColorTexture->transform->texCoordIndex.has_value()) {
					baseColorTexcoordIndex = baseColorTexture->transform->texCoordIndex.value();
				} else {
					baseColorTexcoordIndex = material.pbrData.baseColorTexture->texCoordIndex;
				}
            }
        } else {
			outPrimitive.materialUniformsIndex = 0;
		}
        
            // Position
            auto& positionAccessor = asset.accessors[positionIt->accessorIndex];
            if (!positionAccessor.bufferViewIndex.has_value())
                continue;
			
//			fmt::print("Position Accessor Count: {}\n", positionAccessor.count);
			size_t initVerticesSize = vertices->size();
			vertices->resize(initVerticesSize + positionAccessor.count);
			outPrimitive.vertexStartIdx = initVerticesSize;
			outPrimitive.vertexCount = positionAccessor.count;

			fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec3>(asset, positionAccessor, [&](fastgltf::math::fvec3 pos, std::size_t idx) {
				( *vertices )[ initVerticesSize + idx ].position = glm::vec3(pos.x(), pos.y(), pos.z());
				min = glm::min(min, ( *vertices )[ initVerticesSize + idx ].position);
				max = glm::max(max, ( *vertices )[ initVerticesSize + idx ].position);
				( *vertices )[ initVerticesSize + idx ].normal = glm::vec3(1.0f, 0.0f, 0.0f);
				( *vertices )[ initVerticesSize + idx ].texcoords = glm::vec2(0.0f, 0.0f);
				// fmt::print("Vertex Position: {}\n", glm::to_string(( *vertices )[ initVerticesSize + idx ].position));
			});
        

		auto texcoordAttribute = std::string("TEXCOORD_0");// + std::to_string(baseColorTexcoordIndex);
        if (const auto* texcoord = it->findAttribute(texcoordAttribute); texcoord != it->attributes.end()) {
            // Tex coord
			auto& texCoordAccessor = asset.accessors[texcoord->accessorIndex];
            if (!texCoordAccessor.bufferViewIndex.has_value())
                continue;

			// auto* vertices = static_cast<Vertex*>(glMapNamedBuffer(primitive.vertexBuffer, GL_WRITE_ONLY));
			fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec2>(asset, texCoordAccessor, [&](fastgltf::math::fvec2 uv, std::size_t idx) {
				( *vertices )[ initVerticesSize + idx ].texcoords.x = uv.x();
				( *vertices )[ initVerticesSize + idx ].texcoords.y = uv.y();
			});
        }

		auto normalAttribute = std::string("NORMAL");
        if (const auto* normal = it->findAttribute(normalAttribute); normal != it->attributes.end()) {
            // Tex coord
			auto& normalAccessor = asset.accessors[normal->accessorIndex];
            if (!normalAccessor.bufferViewIndex.has_value())
                continue;

			// auto* vertices = static_cast<Vertex*>(glMapNamedBuffer(primitive.vertexBuffer, GL_WRITE_ONLY));
			fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec3>(asset, normalAccessor, [&](fastgltf::math::fvec3 norm, std::size_t idx) {
				( *vertices )[ initVerticesSize + idx ].normal.x = norm.x();
				( *vertices )[ initVerticesSize + idx ].normal.y = norm.y();
				( *vertices )[ initVerticesSize + idx ].normal.z = norm.z();
			});
		}

		auto& draw = outPrimitive.draw;
		draw.instanceCount = 1;
        draw.baseInstance = 0;
        draw.baseVertex = 0;
		draw.firstIndex = 0;
		
        auto& indexAccessor = asset.accessors[it->indicesAccessor.value()];
		draw.count = static_cast<std::uint32_t>(indexAccessor.count);

		// std::vector<uint> indices (indexAccessor.count, 0u);
		size_t initIndicesSize = indices->size();
		indices->resize( initIndicesSize + indexAccessor.count );
		outPrimitive.indexStartIdx = initIndicesSize;
		outPrimitive.indexCount = indexAccessor.count;


		fastgltf::iterateAccessorWithIndex<uint>(asset, indexAccessor, [&](uint indiceIdx, std::size_t idx) {
				( *indices )[ initIndicesSize + idx ] = indiceIdx;
		});

	
		// Primitive::GenerateBuffers(&outPrimitive, vertices, indices);
		outMesh->primitives.push_back(outPrimitive);
		model->primitives.push_back(outPrimitive);
	    // glVertexArrayElementBuffer(vao, primitive.indexBuffer);

    }

	outMesh->boundingVolume = Sphere((min + max) / 2.0f, glm::distance(min, max) / 2.0f);

    glCreateBuffers(1, &outMesh->drawsBuffer);
    glNamedBufferData(outMesh->drawsBuffer, 
    static_cast<GLsizeiptr>(outMesh->primitives.size() * sizeof(Primitive)),
                      outMesh->primitives.data(), GL_STATIC_DRAW);

	// *outMesh = std::move(retMesh);
    // Create the buffer holding all of our primitive structs.

    return true;
}


bool ModelLoader::LoadImage(const fastgltf::Asset& asset, const fastgltf::Image& image, Texture* outTexture) {	
    auto getLevelCount = [](int width, int height) -> GLsizei {
        return static_cast<GLsizei>(1 + floor(log2(width > height ? width : height)));
    };

    GLuint texture;
    glCreateTextures(GL_TEXTURE_2D, 1, &texture);
    std::visit(fastgltf::visitor {
        [](auto& arg) {
			fmt::print("Loaded Image: {}\n", "Nothing");


		},
        [&](const fastgltf::sources::URI& filePath) {
            assert(filePath.fileByteOffset == 0); // We don't support offsets with stbi.
            assert(filePath.uri.isLocalPath()); // We're only capable of loading local files.
            int width, height, nrChannels;

            const std::string path(filePath.uri.path().begin(), filePath.uri.path().end()); // Thanks C++.
            unsigned char *data = stbi_load(path.c_str(), &width, &height, &nrChannels, 4);
            

			glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
			glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
			glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
			
			glTextureStorage2D(texture, getLevelCount(width, height), GL_RGBA8, width, height);
            glTextureSubImage2D(texture, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);
			fmt::print("Loaded Image: {}\n", path);
        },
        [&](const fastgltf::sources::Array& vector) {
            int width, height, nrChannels;
            unsigned char *data = stbi_load_from_memory(reinterpret_cast<const stbi_uc*>(vector.bytes.data()), static_cast<int>(vector.bytes.size()), &width, &height, &nrChannels, 4);
            
		
			glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
			glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
			glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
			glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

			glTextureStorage2D(texture, getLevelCount(width, height), GL_RGBA8, width, height);
            glTextureSubImage2D(texture, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);
            stbi_image_free(data);
			fmt::print("Loaded Image: {}\n", "From Memory");

        },
        [&](const fastgltf::sources::BufferView& view) {
            auto& bufferView = asset.bufferViews[view.bufferViewIndex];
            auto& buffer = asset.buffers[bufferView.bufferIndex];
            // Yes, we've already loaded every buffer into some GL buffer. However, with GL it's simpler
            // to just copy the buffer data again for the texture. Besides, this is just an example.
            std::visit(fastgltf::visitor {
                // We only care about VectorWithMime here, because we specify LoadExternalBuffers, meaning
                // all buffers are already loaded into a vector.
                [](auto& arg) {},
                [&](const fastgltf::sources::Array& vector) {
                    int width, height, nrChannels;
					unsigned char* data = stbi_load_from_memory(reinterpret_cast<const stbi_uc*>(vector.bytes.data() + bufferView.byteOffset),
					static_cast<int>(bufferView.byteLength), &width, &height, &nrChannels, 4);
                    

					glTextureParameteri(texture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE );
					glTextureParameteri(texture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE );
					glTextureParameteri(texture, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
					glTextureParameteri(texture, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
					
					glTextureStorage2D(texture, getLevelCount(width, height), GL_RGBA8, width, height);
                    glTextureSubImage2D(texture, 0, 0, 0, width, height, GL_RGBA, GL_UNSIGNED_BYTE, data);
                    stbi_image_free(data);
                }
            }, buffer.data);
			fmt::print("Loaded Image: {}\n", "From BufferView");

        },
    }, image.data);

    glGenerateTextureMipmap(texture);
	outTexture->texture = texture;
    return true;

}

bool ModelLoader::LoadMaterial(const fastgltf::Asset& asset, const fastgltf::Material& material, Material* outMaterial) {
    // Material uniforms = {};
    outMaterial->alphaCutoff = material.alphaCutoff;

    outMaterial->baseColorFactor.x = material.pbrData.baseColorFactor.x();
	outMaterial->baseColorFactor.y = material.pbrData.baseColorFactor.y();
	outMaterial->baseColorFactor.z = material.pbrData.baseColorFactor.z();
	outMaterial->baseColorFactor.w = material.pbrData.baseColorFactor.w();
	outMaterial->flags = 0;
    if (material.pbrData.baseColorTexture.has_value()) {
        outMaterial->flags |= MaterialUniformFlags::HasBaseColorTexture;
    }
	fmt::print("Has Base Color Texture: {}\n", outMaterial->flags & MaterialUniformFlags::HasBaseColorTexture);

    // *outMaterial = std::move(uniforms);
    return true;
}

bool ModelLoader::LoadModel(const fastgltf::Asset& asset, Model * model){
	model->textures.resize(asset.images.size());
	model->materials.resize(asset.materials.size() + 1);
	model->meshes.resize(asset.meshes.size());

	for ( size_t imageIdx = 0; imageIdx < asset.images.size(); ++imageIdx ){
		const fastgltf::Image& gltfImage = asset.images[imageIdx];
		Texture& texture = model->textures[imageIdx];
		if (!LoadImage(asset, gltfImage, &texture)){
			assert(false);
			return false;
		}
	}

	// Add a default material
	auto& defaultMaterial = model->materials[0];
	defaultMaterial.baseColorFactor = glm::vec4(1.0f);
	defaultMaterial.alphaCutoff = 0.0f;
	defaultMaterial.flags = 0;

	for ( size_t materialIdx = 1; materialIdx < asset.materials.size() + 1; ++materialIdx ){
		const fastgltf::Material& gltfMaterial = asset.materials[materialIdx - 1];
		Material& material = model->materials[materialIdx];
		if (!LoadMaterial(asset, gltfMaterial, &material)){
			assert(false);
			return false;
		}
		fmt::print("Material: {}\n", materialIdx);
		fmt::print("Alpha Cutoff: {}\n", material.alphaCutoff);
	}

	std::vector<Vertex> vertices;
	std::vector<uint> indices;

	for ( size_t meshIdx = 0; meshIdx < asset.meshes.size(); ++meshIdx ){
		const fastgltf::Mesh& gltfMesh = asset.meshes[meshIdx];
		Mesh& mesh = model->meshes[meshIdx];
		if (!LoadMesh(asset, gltfMesh, model, &mesh, &vertices, &indices)){
			assert(false);
			return false;
		}
	}

	// Create the flat buffer :D
	glCreateBuffers( 1,  &model->VBO );
	glNamedBufferStorage( model->VBO, vertices.size() * sizeof( Vertex ), vertices.data(), GL_DYNAMIC_STORAGE_BIT );

	glCreateBuffers( 1,  &model->EBO );
	glNamedBufferStorage( model->EBO, indices.size() * sizeof( uint ), indices.data(), GL_DYNAMIC_STORAGE_BIT );

	glCreateBuffers(1, &model->nodePropertiesBuffer);
	glNamedBufferStorage( model->nodePropertiesBuffer, model->nodePropertiesBufferVec.size() * sizeof(NodeProperties), model->nodePrimPropertiesBufferVec.data(), GL_DYNAMIC_STORAGE_BIT);
	

	// const auto& sceneIndex = asset.desfaultScene.value_or(0);
	glm::vec3 min = glm::vec3(std::numeric_limits<float>::max());
	glm::vec3 max = glm::vec3(std::numeric_limits<float>::min());
	fastgltf::iterateSceneNodes(asset, 0u, fastgltf::math::fmat4x4(),
								[&](const fastgltf::Node& node, fastgltf::math::fmat4x4 matrix) {
		if (node.meshIndex.has_value()) {
			glm::mat4 modelFromNodeMat = glm::make_mat4(matrix.data());
			model->nodes.push_back( Node( modelFromNodeMat, *node.meshIndex ) );
			min = glm::min(min, model->meshes[*node.meshIndex].boundingVolume.center - glm::vec3(model->meshes[*node.meshIndex].boundingVolume.radius));
			max = glm::max(max, model->meshes[*node.meshIndex].boundingVolume.center + glm::vec3(model->meshes[*node.meshIndex].boundingVolume.radius));

			model->nodePropertiesBufferVec.push_back( { modelFromNodeMat } );
		}
	});

	glCreateBuffers(1, &model->nodePropertiesBuffer);
	glNamedBufferStorage( model->nodePropertiesBuffer, model->nodePropertiesBufferVec.size() * sizeof(NodeProperties), model->nodePropertiesBufferVec.data(), GL_DYNAMIC_STORAGE_BIT );


	model->boundingVolume = Sphere((min + max) / 2.0f, glm::distance(min, max) / 2.0f);
	

	for ( size_t nodeIdx = 0; nodeIdx < model->nodes.size(); ++nodeIdx ) 
	{
		const Mesh& mesh = model->meshes[model->nodes[nodeIdx].meshIndex];
		
		for (size_t primOffset = 0; primOffset < mesh.primCount; ++primOffset)
		{
			const uint primIdx = mesh.startPrimIdx + primOffset;
			model->nodePrimPropertiesBufferVec.push_back( { nodeIdx,  primIdx } );
			const Primitive& prim = model->primitives[ primIdx ];
			model->drawCmdBufferVec.push_back ( {
				prim.indexCount, 		// indexCount
				1,						// instanceCount
				prim.indexStartIdx,		// firstIndex
				prim.vertexStartIdx,	// baseVertex
				0						// baseInstance
			} );
		}
	}
	glCreateBuffers(1, &model->nodePrimPropertiesBuffer);
	glNamedBufferStorage(model->nodePrimPropertiesBuffer, model->nodePrimPropertiesBufferVec.size() * sizeof(NodePrimProperties), model->nodePrimPropertiesBufferVec.data(), GL_DYNAMIC_STORAGE_BIT);
	
	glCreateBuffers(1, &model->drawCmdBuffer);
	glNamedBufferStorage(model->drawCmdBuffer, model->drawCmdBufferVec.size() * sizeof(IndirectDrawCommand), model->drawCmdBufferVec.data(), GL_DYNAMIC_STORAGE_BIT);
	
	
	
	fmt::print("AlphaCutoff after model load {}\n", model->materials[0].alphaCutoff);
	return true;
}

// void ModelLoader::PrintModel(const Model& model){

// }
