#version 460 core

out vec4 FragColor;


layout (location = 3) in vec2 TexCoord;


uniform sampler2D albedoTexture;

void main()
{
    FragColor = texture(albedoTexture, TexCoord);
}
