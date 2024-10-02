#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;
// layout (location = 3) in mat4 worldFromModel;

// layout (location = 3) in int aSelected;

layout (binding = 0, std140) uniform Matrices
{
    mat4 projection;
    mat4 view;    
};

layout (binding = 1, std140) uniform meshPropertiesBuf {
    mat4 modelFromMesh;
};

// SSBO containing the instanced model matrices
layout(binding = 2, std430) readonly buffer instModelMatBuf {
    mat4 worldFromModel[];
};

layout(binding = 3, std430) coherent buffer instMeshRenderedBuf {
    bool meshIsRendered[];
}; // THIS DOES NOT SAVE COMPUTE!!!


layout (location = 10) out vec3 Normal;
layout (location = 11) out vec2 TexCoord;

void main()
{
    if ( true ) // meshIsRendered[gl_InstanceID] )
    {
        gl_Position = projection * view * worldFromModel[gl_DrawID] * modelFromMesh * vec4(aPos, 1.0);
        Normal = aNormal;
        TexCoord =aTexCoord;
    } else {
        gl_Position = vec4(0.0, 0.0, 0.0, 0.0);
        Normal = vec3(0.0, 0.0, 0.0);
        TexCoord = vec2(0.0, 0.0);
    }

}
