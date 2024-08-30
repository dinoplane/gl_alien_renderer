#version 460 core

out vec4 FragColor;

in vec3 Normal;
in vec2 TexCoord;

void main()
{
    vec3 v_normalColor = (normalize(Normal) * 0.25f ) + 0.25f;
    // FragColor = vec4(TexCoord*0.25f, 0.0, 1.0) + vec4(v_normalColor, 1.0);
    // FragColor = vec4(TexCoord*0.5f, 0.0, 1.0);
    FragColor = vec4(v_normalColor, 1.0);
}
