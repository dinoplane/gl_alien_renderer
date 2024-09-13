 #include <model_loader.hpp>

 #include <iostream>
 #include <util.h>

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
			fastgltf::Extensions::KHR_materials_variants;

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

bool ModelLoader::LoadMesh(const fastgltf::Asset& asset, const fastgltf::Mesh& mesh, 
Mesh* outMesh) {
    // Mesh outMesh = {};
	std::vector<Vertex> vertices;
	vertices.resize()
    for (auto it = mesh.primitives.begin(); it != mesh.primitives.end(); ++it) {
		auto* positionIt = it->findAttribute("POSITION");
		assert(positionIt != it->attributes.end()); // A mesh primitive is required to hold the POSITION attribute.
		assert(it->indicesAccessor.has_value()); // We specify GenerateMeshIndices, so we should always have indices

        // // Generate the VAO
        // GLuint vao = GL_NONE;
        // glCreateVertexArrays(1, &vao);

		std::size_t baseColorTexcoordIndex = 0;

        // Get the output primitive
        auto index = std::distance(mesh.primitives.begin(), it);
        // auto& primitive = outMesh.primitives[index];
        // primitive.primitiveType = fastgltf::to_underlying(it->type);
        // primitive.vertexArray = vao;
        // if (it->materialIndex.has_value()) {
        //     primitive.materialUniformsIndex = it->materialIndex.value() + 1; // Adjust for default material
        //     auto& material = viewer->asset.materials[it->materialIndex.value()];

		// 	auto& baseColorTexture = material.pbrData.baseColorTexture;
        //     if (baseColorTexture.has_value()) {
        //         auto& texture = viewer->asset.textures[baseColorTexture->textureIndex];
		// 		if (!texture.imageIndex.has_value())
		// 			return false;
        //         primitive.albedoTexture = viewer->textures[texture.imageIndex.value()].texture;

		// 		if (baseColorTexture->transform && baseColorTexture->transform->texCoordIndex.has_value()) {
		// 			baseColorTexcoordIndex = baseColorTexture->transform->texCoordIndex.value();
		// 		} else {
		// 			baseColorTexcoordIndex = material.pbrData.baseColorTexture->texCoordIndex;
		// 		}
        //     }
        // } else {
			// primitive.materialUniformsIndex = 0;
		// }

        {
            // Position
            auto& positionAccessor = asset.accessors[positionIt->accessorIndex];
            if (!positionAccessor.bufferViewIndex.has_value())
                continue;

			// Create the vertex buffer for this primitive, and use the accessor tools to copy directly into the mapped buffer.
			// glCreateBuffers(1, &primitive.vertexBuffer);
			// glNamedBufferData(primitive.vertexBuffer, positionAccessor.count * sizeof(Vertex), nullptr, GL_STATIC_DRAW);
			// auto* vertices = static_cast<Vertex*>(glMapNamedBuffer(primitive.vertexBuffer, GL_WRITE_ONLY));
			fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec3>(asset, positionAccessor, [&](fastgltf::math::fvec3 pos, std::size_t idx) {
				vertices[idx].position = glm::vec3(pos.x(), pos.y(), pos.z());
				vertices[idx].normal = glm::vec3(0.0f, 0.0f, 0.0f);
				vertices[idx].texcoords = glm::vec2(0.0f, 0.0f);
			});
			// glUnmapNamedBuffer(primitive.vertexBuffer);

            // glEnableVertexArrayAttrib(vao, 0);
            // glVertexArrayAttribFormat(vao, 0,
            //                           3, GL_FLOAT,
            //                           GL_FALSE, 0);
            // glVertexArrayAttribBinding(vao, 0, 0);

			// glVertexArrayVertexBuffer(vao, 0, primitive.vertexBuffer,
			// 						  0, sizeof(Vertex));
        }

		auto texcoordAttribute = std::string("TEXCOORD_") + std::to_string(baseColorTexcoordIndex);
        if (const auto* texcoord = it->findAttribute(texcoordAttribute); texcoord != it->attributes.end()) {
            // Tex coord
			auto& texCoordAccessor = asset.accessors[texcoord->accessorIndex];
            if (!texCoordAccessor.bufferViewIndex.has_value())
                continue;

			// auto* vertices = static_cast<Vertex*>(glMapNamedBuffer(primitive.vertexBuffer, GL_WRITE_ONLY));
			fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec2>(asset, texCoordAccessor, [&](fastgltf::math::fvec2 uv, std::size_t idx) {
				vertices[idx].texcoords.x = uv.x();
				vertices[idx].texcoords.y = uv.y();
			});
			// glUnmapNamedBuffer(primitive.vertexBuffer);

			// glEnableVertexArrayAttrib(vao, 1);
            // glVertexArrayAttribFormat(vao, 1,
			// 						  2, GL_FLOAT,
            //                           GL_FALSE, 0);
            // glVertexArrayAttribBinding(vao, 1, 1);

			// glVertexArrayVertexBuffer(vao, 1, primitive.vertexBuffer,
			// 						  offsetof(Vertex, uv), sizeof(Vertex));
        }

		auto normalAttribute = std::string("NORMAL");
        if (const auto* normal = it->findAttribute(normalAttribute); normal != it->attributes.end()) {
            // Tex coord
			auto& normalAccessor = asset.accessors[normal->accessorIndex];
            if (!normalAccessor.bufferViewIndex.has_value())
                continue;

			// auto* vertices = static_cast<Vertex*>(glMapNamedBuffer(primitive.vertexBuffer, GL_WRITE_ONLY));
			fastgltf::iterateAccessorWithIndex<fastgltf::math::fvec3>(asset, normalAccessor, [&](fastgltf::math::fvec3 norm, std::size_t idx) {
				vertices[idx].normal.x = norm.x();
				vertices[idx].normal.y = norm.y();
				vertices[idx].normal.z = norm.z();
			});
			// glUnmapNamedBuffer(primitive.vertexBuffer);

			// glEnableVertexArrayAttrib(vao, 1);
            // glVertexArrayAttribFormat(vao, 1,
			// 						  2, GL_FLOAT,
            //                           GL_FALSE, 0);
            // glVertexArrayAttribBinding(vao, 1, 1);

			// glVertexArrayVertexBuffer(vao, 1, primitive.vertexBuffer,
			// 						  offsetof(Vertex, uv), sizeof(Vertex));
        }

        // Generate the indirect draw command
        // auto& draw = primitive.draw;
        // draw.instanceCount = 1;
        // draw.baseInstance = 0;
        // draw.baseVertex = 0;
		// draw.firstIndex = 0;


        auto& indexAccessor = asset.accessors[it->indicesAccessor.value()];
        if (!indexAccessor.bufferViewIndex.has_value())
            return false;
        // draw.count = static_cast<std::uint32_t>(indexAccessor.count);
		

		// // Create the index buffer and copy the indices into it.
		// // glCreateBuffers(1, &primitive.indexBuffer);
		// if (indexAccessor.componentType == fastgltf::ComponentType::UnsignedByte || indexAccessor.componentType == fastgltf::ComponentType::UnsignedShort) {
        // 	// primitive.indexType = GL_UNSIGNED_SHORT;
		// 	// glNamedBufferData(primitive.indexBuffer,
		// 	// 				  static_cast<GLsizeiptr>(indexAccessor.count * sizeof(std::uint16_t)), nullptr,
		// 	// 				  GL_STATIC_DRAW);
		// 	// auto* indices = static_cast<std::uint16_t*>(glMapNamedBuffer(primitive.indexBuffer, GL_WRITE_ONLY));
		// 	uint16_t* indexPtr;
		// 	fastgltf::copyFromAccessor<std::uint16_t>(asset, indexAccessor, indexPtr);
		// 	// glUnmapNamedBuffer(primitive.indexBuffer);
		// } else {
        // 	// primitive.indexType = GL_UNSIGNED_INT;
		// 	// glNamedBufferData(primitive.indexBuffer,
		// 					//   static_cast<GLsizeiptr>(indexAccessor.count * sizeof(std::uint32_t)), nullptr,
		// 					//   GL_STATIC_DRAW);
		// 	// auto* indices = static_cast<std::uint32_t*>(glMapNamedBuffer(primitive.indexBuffer, GL_WRITE_ONLY));
		// 	uint32_t* indexPtr;
		// 	fastgltf::copyFromAccessor<std::uint32_t>(asset, indexAccessor, indexPtr);
		// 	// glUnmapNamedBuffer(primitive.indexBuffer);
		// }


		uint* indexPtr;
		fastgltf::copyFromAccessor<uint>(asset, indexAccessor, indexPtr);
		std::vector<uint> indices (indexPtr, indexPtr + indexAccessor.count);

	Mesh::GenerateBuffers(outMesh, vertices, indices);
        // glVertexArrayElementBuffer(vao, primitive.indexBuffer);
    }

    // Create the buffer holding all of our primitive structs.

    return true;
}

/*
void drawMesh(Viewer* viewer, std::size_t meshIndex, fastgltf::math::fmat4x4 matrix) {
    auto& mesh = viewer->meshes[meshIndex];

    glBindBuffer(GL_DRAW_INDIRECT_BUFFER, mesh.drawsBuffer);

    glUniformMatrix4fv(viewer->modelMatrixUniform, 1, GL_FALSE, &matrix[0][0]);

    for (auto i = 0U; i < mesh.primitives.size(); ++i) {
        auto& prim = mesh.primitives[i];
		auto& gltfPrimitive = viewer->asset.meshes[meshIndex].primitives[i];

		std::size_t materialIndex;
		auto& mappings = gltfPrimitive.mappings;
		if (!mappings.empty() && mappings[viewer->materialVariant].has_value()) {
			materialIndex = mappings[viewer->materialVariant].value() + 1; // Adjust for default material
		} else {
			materialIndex = prim.materialUniformsIndex;
		}

        auto& material = viewer->materialBuffers[materialIndex];
        glBindTextureUnit(0, prim.albedoTexture);
        glBindBufferBase(GL_UNIFORM_BUFFER, 0, material);
        glBindVertexArray(prim.vertexArray);

		// Update texture transform uniforms
		`glUniform2f(viewer->uvOffsetUniform, 0, 0);
		glUniform2f(viewer->uvScaleUniform, 1.f, 1.f);
		glUniform1f(viewer->uvRotationUniform, 0);
		if (materialIndex != 0) {
			auto& gltfMaterial = viewer->asset.materials[materialIndex - 1];
			if (gltfMaterial.pbrData.baseColorTexture.has_value() && gltfMaterial.pbrData.baseColorTexture->transform) {
				auto& transform = gltfMaterial.pbrData.baseColorTexture->transform;
				glUniform2f(viewer->uvOffsetUniform, transform->uvOffset[0], transform->uvOffset[1]);
				glUniform2f(viewer->uvScaleUniform, transform->uvScale[0], transform->uvScale[1]);
				glUniform1f(viewer->uvRotationUniform, static_cast<float>(transform->rotation));
			}
		}

        glDrawElementsIndirect(prim.primitiveType, prim.indexType,
                               reinterpret_cast<const void*>(i * sizeof(Primitive)));
    }
}

	// Load the glTF file
    auto start = std::chrono::high_resolution_clock::now();
    if (!loadGltf(&viewer, gltfFile)) {
        std::cerr << "Failed to parse glTF" << '\n';
        return -1;
    }

    for (auto& mesh : asset.meshes) {
        loadMesh(&viewer, mesh);
    }


*/