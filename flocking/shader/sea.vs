#version 330 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;
// layout (location = 3) in int aSelected;

#define WAVECOUNT 4

#define PI 3.1415926535897932384626433832795


uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform int selected;
uniform float uTime;


out vec2 TexCoord;
flat out int Selected;
out vec4 FragPos;
out vec3 Normal;


struct Wave{
    float L; // Wavelength, frequency w = 2pi/L
    float A; // Amplitude
    float S; // Speed, phase constant phi = S * 2pi/L
    vec2 D; // Direction of the wave with respect to 2d surface
};


Wave waveArray[WAVECOUNT] = Wave[WAVECOUNT](
    Wave(5.0, 0.1, 5, vec2(1.0, 0.0)),
    Wave(7.0, 0.05, 1, vec2(-1.0, 1.0)),
    Wave(10.0, 0.1, 4, vec2(0.0, 1.0)),
    Wave(8.0, 0.3, 2, vec2(-1.0, -1.0))
);

void main()
{
    // gl_Position = view  * vec4(aPos, 1.0);

    vec4 pos = model * vec4(aPos, 1.0);
    // if (pos.x > 0.0)
    //     pos.y += 0.5;
    float x = pos.x;
    float z = pos.z;
    float t = uTime;

    Normal = vec3(0, 1, 0);
    // General sum of sines
    for (int i = 0; i < WAVECOUNT; i++){
        Wave wv = waveArray[i];
        float w = 2*PI / wv.L;
        pos.y += wv.A * sin( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) );
        Normal.x -= w * wv.D.x * wv.A * cos( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) );
        Normal.z -= w * wv.D.y * wv.A * cos( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) );
    }

    // pos.y += sin(0.5*i) + sin(i);

    // Normal = vec3(0.5*cos(0.5*i), 1.0, 0.0);
    // vec3 t = vec3(1.0, 0.5*cos(0.5*i) + cos(i), 0.0);
    // vec3 b = vec3(0.0, 0.0, 1.0);

    // Normal = cross(b, t);

    FragPos = pos;
    gl_Position = projection * view * FragPos;


    Selected = selected;
    // ourColor = aPos + vec3(0.5);
    TexCoord = vec2(aTexCoord.x, 1.0 - aTexCoord.y);

}
