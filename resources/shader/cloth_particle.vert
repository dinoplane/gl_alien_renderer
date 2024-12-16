#version 460 core
layout (location = 0) in vec4 aPos;


layout (location = 3) out vec2 TexCoord;

layout (std140, binding=0) uniform Matrices
{
    mat4 projection;
    mat4 view;
};



// SSBO containing the texcoords of particles
layout(binding = 8, std430) buffer texcoordBuffer {
    vec2 texcoords[]; // arrlen total # of particles, indexed with gl_VertexID
};

void main()
{
    gl_Position = projection * view * aPos;
    TexCoord = texcoords[gl_VertexID];
}
