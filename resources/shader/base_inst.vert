#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aNormal;
layout (location = 2) in vec2 aTexCoord;
layout (location = 3) in mat4 worldFromModel;


// layout (location = 3) in int aSelected;

uniform mat4 view;
uniform mat4 projection;

out vec3 Normal;
out vec2 TexCoord;

void main()
{
    gl_Position = projection * view * worldFromModel * vec4(aPos, 1.0);
    Normal = aNormal;
    TexCoord = vec2(aPos.x, 1.0 - aPos.y);
}