#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;
// layout (location = 3) in mat4 worldFromModel;

// layout (location = 3) in int aSelected;

layout (std140, binding=3) uniform Matrices
{
    mat4 projection;
    mat4 view;
};


// SSBO containing the instanced model matrices
layout(binding = 4, std430) readonly buffer instModelMatBuf {
    mat4 worldFromModel[];
};

layout(binding = 5, std430) coherent buffer instMeshRenderedBuf {
    bool meshIsRendered[];
}; // THIS DOES NOT SAVE COMPUTE!!!

layout (binding = 6, std140) uniform meshPropertiesBuf {
    mat4 modelFromMesh;
};


out vec3 Normal;
out vec2 TexCoord;

void main()
{
    if ( meshIsRendered[gl_InstanceID] )
    {
        gl_Position = projection * view * worldFromModel[gl_InstanceID] * modelFromMesh * vec4(aPos, 1.0);
        Normal = aNormal;
        TexCoord = vec2(aPos.x, 1.0 - aPos.y);
    } else {
        gl_Position = vec4(0.0, 0.0, 0.0, 0.0);
        Normal = vec3(0.0, 0.0, 0.0);
        TexCoord = vec2(0.0, 0.0);
    }

}
