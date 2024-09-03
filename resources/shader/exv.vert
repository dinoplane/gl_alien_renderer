#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;
// layout (location = 3) in int aSelected;


uniform mat4 model;
layout (std430, binding=3) uniform Matrices
{
    mat4 projection;
    mat4 view;
};
uniform int selected;


out vec2 TexCoord;
flat out int Selected;

void main()
{
    // gl_Position = view  * vec4(aPos, 1.0);

    gl_Position = projection * view * model * vec4(aPos, 1.0);


    Selected = selected;
    // ourColor = aPos + vec3(0.5);
    TexCoord = vec2(aTexCoord.x, 1.0 - aTexCoord.y);
}
