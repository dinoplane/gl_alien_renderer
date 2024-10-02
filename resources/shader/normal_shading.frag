#version 460 core

out vec4 FragColor;

layout (location = 10) in vec3 Normal;
layout (location = 11) in vec2 TexCoord;

void main()
{
    vec3 v_normalColor = (normalize(Normal) * 0.25f ) + 0.25f;
    FragColor = vec4(v_normalColor, 1.0);
}
