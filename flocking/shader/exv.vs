#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;

out vec3 ourColor;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

out vec2 TexCoord;

void main()
{
    // gl_Position = view  * vec4(aPos, 1.0);

    gl_Position = projection * view * model * vec4(aPos, 1.0);

    if (aPos.x == 0.5){
        ourColor = vec3(1.0, 0.0, 0.0);
    } else if (aPos.y == 0.5){
        ourColor = vec3(0.0, 1.0, 0.0);
    } else if (aPos.z == 0.5){
        ourColor = vec3(0.0, 0.0, 1.0);
    } else if (aPos.x == -0.5){
        ourColor = vec3(0.0, 1.0, 1.0);
    } else if (aPos.y == -0.5){
        ourColor = vec3(1.0, 0.0, 1.0);
    } else if (aPos.z == -0.5){
        ourColor = vec3(1.0, 1.0, 0.0);
    }
    // ourColor = aPos + vec3(0.5);
    TexCoord = vec2(aTexCoord.x, 1.0 - aTexCoord.y);
}
