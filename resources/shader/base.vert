#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;

uniform mat4 model;
layout (std140, binding=0) uniform Matrices
{
    mat4 projection;
    mat4 view;
};

layout (location = 10) out vec3 Normal;
layout (location = 11) out vec2 TexCoord;

void main()
{
    gl_Position = projection * view * model * vec4(aPos, 1.0);
    Normal = aNormal;
    TexCoord = vec2(aPos.x, 1.0 - aPos.y);
}
