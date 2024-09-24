#version 460 core
out vec4 FragColor;

layout (location = 10) in vec2 TexCoords;
uniform sampler2D screenTexture;

void main()
{
    FragColor = texture(screenTexture, TexCoords.st);
}