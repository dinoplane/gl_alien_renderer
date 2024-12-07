#version 460 core

out vec4 FragColor;

layout (location = 10) in vec3 Normal;
layout (location = 11) in vec2 TexCoord;
layout (location = 12) in flat uint materialIdx;

uniform sampler2D albedoTexture;

struct Material {
    vec4 baseColorFactor;
    float alphaCutoff;
    uint flags;
};

float rand(vec2 co){
    return fract(sin(dot(co, vec2(12.9898, 78.233))) * 43758.5453);
}


layout (binding = 6, std430) readonly buffer materialPropertiesBuf {
    Material materialProperties[]; // arrlen total materials in model, indexed with materialIndices[gl_DrawID]
};

void main()
{
    vec4 color = texture(albedoTexture, TexCoord);

    Material material = materialProperties[materialIdx];

    float factor = (rand(gl_FragCoord.xy) - 0.5) / 8;
    if (color.a < material.alphaCutoff)
        discard;
    else {
        // color.a = 1.0;
    }

    FragColor = color;
}
