#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;

// Structs -------------------------------------------------------------
struct NodePrimProperties {
    uint nodeIdx;
    uint primIdx;
};

struct NodeProperties {
    mat4 modelFromNode;
};

struct Material {
    vec4 baseColorFactor;
    float alphaCutoff;
    uint flags;
};

struct Primitive {
    uint materialIdx;
    uint textureIdx;
};

struct AtlasOffsetsData {
    float offsetX;
    float mipLevel; // Unused, but there for padding
    float textureWidth;
    float textureHeight;
};

// ---------------------------------------------------------------------



// Uniforms ------------------------------------------------------------

// UBO containing the projection and view matrices
layout (binding = 0, std140) uniform Matrices
{
    mat4 projection;
    mat4 view;
};

// ---------------------------------------------------------------------


// Shader Storage Buffers ----------------------------------------------

// SSBO containing a mapping between nodes and primitives
layout (binding = 1, std430) readonly buffer nodePrimPropertiesBuf {
    NodePrimProperties nodePrimProperties[]; // arrlen total logical primitives in model, indexed with gl_DrawID
};

// SSBO containing the primitive properties
layout (binding = 2, std430) readonly buffer primitivePropertiesBuf {
    Primitive primitiveProperties[]; // arrlen total primitives in model, indexed with nodePrimProperties[gl_DrawID].primIdx
};


// SSBO containing the indices of the visible instances
layout(binding = 3, std430) readonly buffer visibleInstIndicesBuf {
    uint visibleInstCount;
    uint visibleInstIndices[]; // arrlen total visible instances in scene; indexed with gl_InstanceID
};

// SSBO containing the instanced model matrices
layout(binding = 4, std430) readonly buffer instModelMatBuf {
    mat4 worldFromModel[]; // arrlen total instances in scene; indexed with visibleInstIndices[gl_InstanceID]
};

// SSBO containing the node matrices
layout(binding = 5, std430) readonly buffer NodePropertiesBuf {
    NodeProperties nodeProperties[]; // arrlen total nodes in model, indexed with nodePrimProperties[gl_DrawID].nodeIdx
};


// Note: The below are only relevant in the fragment shader, pass the indices to the fragment shader

// SSBO containing the material properties
layout (binding = 6, std430) readonly buffer materialPropertiesBuf {
    Material materialProperties[]; // arrlen total materials in model, indexed with materialIndices[gl_DrawID]
};


layout (binding = 7, std430) readonly buffer atlasOffsetsBuf {
    float toffsetX;
    float tmipLevel; // Unused, but there for padding
    float totalTextureWidth;
    float totalTextureHeight;
    AtlasOffsetsData atlasOffsets[]; // arrlen total materials in model, indexed with materialIndices[gl_DrawID]
};


// SSBO containing the texture properties
// layout (binding = 7, std430) readonly buffer materialPropertiesBuf {
//     sampler2D albedoTextures[]; // arrlen total textures in model, indexed with textureIndices[gl_DrawID]
// };

// ---------------------------------------------------------------------



/*
struct inputData {
    mat4 worldFromModel; // = worldFromModel[visibleInstIndices[gl_InstanceID]];
    mat4 modelFromMesh; // = nodeProperties[nodePrimProperties[gl_DrawID].nodeIdx]
    Material material; // = materialProperties[primitiveProperties[nodePrimProperties[gl_DrawID].primIdx].materialIdx];
    sampler2D tex; // = primitiveProperties[nodePrimProperties[gl_DrawID].primIdx].textureIdx;
};

*/

layout (location = 10) out vec3 Normal;
layout (location = 11) out vec2 TexCoord;
layout (location = 12) out uint materialIdx;

void main()
{
    gl_Position = projection * view * worldFromModel[visibleInstIndices[gl_InstanceID]] * nodeProperties[nodePrimProperties[gl_DrawID].nodeIdx].modelFromNode * vec4(aPos, 1.0);
    Normal = aNormal;
    float offsetX = atlasOffsets[primitiveProperties[nodePrimProperties[gl_DrawID].primIdx].textureIdx].offsetX;
    float textureWidth = atlasOffsets[primitiveProperties[nodePrimProperties[gl_DrawID].primIdx].textureIdx].textureWidth;
    float textureHeight = atlasOffsets[primitiveProperties[nodePrimProperties[gl_DrawID].primIdx].textureIdx].textureHeight;

    TexCoord.x = (offsetX + aTexCoord.x * textureWidth) / totalTextureWidth;
    TexCoord.y = aTexCoord.y * textureHeight / totalTextureHeight;
    materialIdx = primitiveProperties[nodePrimProperties[gl_DrawID].primIdx].materialIdx;
}
