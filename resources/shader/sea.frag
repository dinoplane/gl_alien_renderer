#version 330 core
out vec4 FragColor;

flat in int Selected;
in vec2 TexCoord;
in vec3 FragPos;
in vec3 Normal;


uniform sampler2D texture1;
uniform sampler2D texture2;
uniform vec3 uLightPos;
uniform vec3 uViewPos;

void main()
{
    // Diffuse
    // vec4 diffColor = vec4(0.1, 0.2, 0.4, 1.0);
    vec4 diffColor = vec4(0.62, 0.85, 0.84, 1.0);

    vec3 lightVec = normalize(uLightPos - FragPos);
    vec3 norm = normalize(Normal);
    float diff = max(dot(lightVec, norm), 0.0);
    if (diff  == 0 && uLightPos.y <= 0){
        diff = max(dot(lightVec, -norm), 0.0);
    }
    float albedo = 1.0;
    vec4 diffuse = albedo * diff * diffColor;

    // Specular
    float k = 32;
    float sp = 0.5;
    vec4 specColor = vec4(1.0);
    vec3 viewDir = normalize(uViewPos - FragPos);
    vec3 reflectDir = reflect(-lightVec, norm);
    vec4 specular = specColor * sp * pow(max(dot(viewDir, reflectDir), 0), k);

    // Ambient
    float ambience = 0.0;
    vec4 ambColor = vec4(1.0, 0.0, 0.0, 1.0);
    vec4 ambient = ambColor * ambience;

    // FragColor = vec4(ourColor, 1.0);
    if (Selected == 0){
        FragColor = diffuse + specular + ambient;
    } else FragColor = mix(texture(texture1, TexCoord), texture(texture2, TexCoord), 0.2);
}
