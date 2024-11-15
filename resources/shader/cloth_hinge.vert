#version 460 core
layout (location = 0) in vec4 aPos;

layout (std140, binding=0) uniform Matrices
{
    mat4 projection;
    mat4 view;
};


layout (std140, binding=1) uniform pickedDebugData
{
    uint pickedEdgeIdx;
    uint pickedHingeIdx;
};


// SSBO containing the positions of the particles
layout(binding = 0, std430) buffer particlePositionsBuf {
    vec4 positions[]; // arrlen total # of particles, indexed with globalInvocationID.x
};

// SSBO containing the edges of the cloth, which are indices into the particle positions
layout(binding = 6, std430) buffer particleHingeBuf {
    ivec4 hinges[]; // arrlen total # of hinges, indexed with pickedHingeIdx
                    // each edge indexed with gl_VertexID
};

void main()
{
    gl_Position = projection * view * positions[hinges[pickedHingeIdx][gl_VertexID]];
}
