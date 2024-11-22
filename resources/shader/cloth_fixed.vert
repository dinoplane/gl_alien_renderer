#version 460 core
layout (location = 0) in uint aIdx;

layout (std140, binding=0) uniform Matrices
{
    mat4 projection;
    mat4 view;
};

// SSBO containing the positions of the particles
layout(binding = 0, std430) buffer particlePositionsBuf {
    vec4 positions[]; // arrlen total # of particles, indexed with globalInvocationID.x
};



void main()
{
    gl_Position = projection * view * positions[aIdx];
}
