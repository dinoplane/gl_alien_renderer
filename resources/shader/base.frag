#version 460 core

out vec4 FragColor;

layout (location = 0) in vec3 Normal;
layout (location = 1) in vec2 TexCoord;


const uint HAS_BASE_COLOR_TEXTURE = 1;

layout(location = 4) uniform sampler2D albedoTexture;
layout(binding = 5, std140) uniform MaterialUniformsUBO {
    vec4 baseColorFactor;
    float alphaCutoff;
    uint flags;
} material;

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}

void main()
{
    // vec3 v_normalColor = (normalize(Normal) * 0.25f ) + 0.25f;
    // // FragColor = vec4(TexCoord*0.25f, 0.0, 1.0) + vec4(v_normalColor, 1.0);
    // // FragColor = vec4(TexCoord*0.5f, 0.0, 1.0);
    // FragColor = vec4(v_normalColor, 1.0);

    vec4 color = material.baseColorFactor;
    if ((material.flags & HAS_BASE_COLOR_TEXTURE) == HAS_BASE_COLOR_TEXTURE) {
        color *= texture(albedoTexture, TexCoord);

        // color *= texture(albedoTexture, transformUv(texCoord));
    }
    float factor = (rand(gl_FragCoord.xy) - 0.5) / 8;
    if (color.a < material.alphaCutoff + factor)
        discard;

    FragColor = color;

}
