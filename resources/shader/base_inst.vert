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
layout(binding = 4, std430) readonly buffer ssbo1 {
    mat4 worldFromModel[];
};


out vec3 Normal;
out vec2 TexCoord;

void main()
{
    gl_Position = projection * view * worldFromModel[gl_InstanceID] * vec4(aPos, 1.0);
    Normal = aNormal;
    TexCoord = vec2(aPos.x, 1.0 - aPos.y);
}
