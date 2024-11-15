#version 460 core
layout (location = 0) in vec4 aPos;

layout (std140, binding=0) uniform Matrices
{
    mat4 projection;
    mat4 view;
};

void main()
{
    gl_Position = projection * view * aPos;
}
