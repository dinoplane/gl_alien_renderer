#version 330 core
out vec4 FragColor;

flat in int Selected;
in vec2 TexCoord;
in vec4 FragPos;
in vec3 Normal;


uniform sampler2D texture1;
uniform sampler2D texture2;
uniform vec3 uLightPos;

void main()
{

    vec3 lightVec = normalize(uLightPos - vec3(FragPos));
    vec3 norm = normalize(Normal);
    float diff = max(dot(lightVec, norm), 0.0);
    float albedo = 1.0;

    vec4 color = vec4(1.0, 0.0, 0.0, 1.0);
    // FragColor = vec4(ourColor, 1.0);
    if (Selected == 0){
        FragColor = albedo*diff*color;
    } else FragColor = mix(texture(texture1, TexCoord), texture(texture2, TexCoord), 0.2);
}
