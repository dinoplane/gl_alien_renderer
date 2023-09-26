#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;
// layout (location = 3) in int aSelected;


uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform int selected;
uniform float uTime;


out vec2 TexCoord;
flat out int Selected;
out vec4 FragPos;
out vec3 Normal;

void main()
{
    // gl_Position = view  * vec4(aPos, 1.0);

    vec4 pos = model * vec4(aPos, 1.0);
    // if (pos.x > 0.0)
    //     pos.y += 0.5;
    float i = pos.x + uTime;
    float j = pos.z + uTime;
    pos.y += sin(0.5*i) + sin(i) + sin(0.3*j);

    // Normal = vec3(0.5*cos(0.5*i), 1.0, 0.0);
    vec3 t = vec3(1.0, 0.5*cos(0.5*i) + cos(i), 0.0);
    vec3 b = vec3(0.0, 0.3*cos(0.3*j), 1.0);

    Normal = cross(b, t);

    FragPos = pos;
    gl_Position = projection * view * FragPos;


    Selected = selected;
    // ourColor = aPos + vec3(0.5);
    TexCoord = vec2(aTexCoord.x, 1.0 - aTexCoord.y);

}
