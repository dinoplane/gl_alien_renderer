#version 460 core
layout (triangles) in;
layout (line_strip, max_vertices = 6) out;

in VS_OUT {
    vec4 normal;
} gs_in[];

const float MAGNITUDE = 0.4;
layout (std430, binding=3) uniform Matrices
{
    mat4 projection;
    mat4 view;
};

void GenerateLine(int index)
{
    gl_Position =  gl_in[index].gl_Position;
    EmitVertex();
    gl_Position = (gl_in[index].gl_Position +
                gs_in[index].normal * MAGNITUDE);
    EmitVertex();
    EndPrimitive();
}

void main()
{
    GenerateLine(0); // first vertex normal
    GenerateLine(1); // second vertex normal
    GenerateLine(2); // third vertex normal
}