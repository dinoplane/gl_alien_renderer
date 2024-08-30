#version 460 core
out vec4 FragColor;

flat in int Selected;
in vec2 TexCoord;

in vec3 Position;

uniform sampler2D noise1;

uniform sampler2D texture1;
uniform sampler2D texture2;

float noise (in vec2 p){
    vec4 val = texture(noise1, p);
    return val.x;
}

void main()
{
    // FragColor = vec4(ourColor, 1.0);
    if (Selected == 0){
        FragColor = vec4(1.0, 0.0, 0.0, 1.0);
    } else FragColor = vec4(((Position) + 1.0)/2.0, 1.0);
    // if (Selected == 1){
    //     FragColor = vec4(1.0, 0.0, 0.0, 1.0);
    // } else FragColor = mix(texture(texture1, TexCoord), texture(texture2, TexCoord), 0.2);
}
