#version 460 core
layout (location = 0) in vec3 aPos;
layout (location = 1) in vec3 aColor;
layout (location = 2) in vec2 aTexCoord;
// layout (location = 3) in int aSelected;

#define WAVECOUNT 4

#define PI 3.1415926535897932384626433832795

#define K_C 4.5

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;
uniform int selected;
uniform float uTime;

out vec2 TexCoord;
flat out int Selected;
out vec3 FragPos;
out vec3 Normal;
out vec3 oldPos;
out vec3 newPos;


float random (vec2 st) {
    return fract(sin(dot(st.xy,
                         vec2(12.9898,78.233)))*
        43758.5453123);
}

float random (float i) {
    return random(vec2(i));
}

struct Wave{
    float L; // Wavelength, frequency w = 2pi/L
    float A; // Amplitude
    float S; // Speed, phase constant phi = S * 2pi/L
    vec2 D; // Direction of the wave with respect to 2d surface
};


Wave waveArray[WAVECOUNT] = Wave[WAVECOUNT](
    Wave(14.0, 0.07, 5, vec2(1.0, 5.0)),
    Wave(7.0, 0.05, 1, vec2(-1.0, 1.0)),
    Wave(10.0, 0.07, 4, vec2(1.0, 0.5)),
    Wave(8.0, 0.15, 2, vec2(-1.0, -1.0))
);

float waveFunction(float x, float z, float t, Wave wv){
    float w = 2*PI / wv.L;
    return wv.A * sin( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) );
    // return 2.0 * wv.A * pow(((sin( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) ) + 1) / 2), K_C) ;
}

float waveDerivativeX(float x, float z, float t, Wave wv){
    float w = 2*PI / wv.L;
    return -w * wv.D.x * wv.A * cos( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) );
    // return -K_C * w * wv.D.x * wv.A * cos( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) )
    //         * pow(((sin( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) ) + 1) / 2), (K_C - 1));

}
float waveDerivativeZ(float x, float z, float t, Wave wv){
    float w = 2*PI / wv.L;
    return -w * wv.D.y * wv.A * cos( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) );

    // return -K_C * w * wv.D.y * wv.A * cos( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) )
    //         * pow(((sin( dot(wv.D, vec2(x, z) * w) + t * (wv.S * w) ) + 1) / 2), (K_C - 1));

}


// GPU Gems 1 Ch. 2
vec3 line_plane_intercept(vec3 POI, vec3 planeP, vec3 planeN, vec3 lineU) {
    float d = dot(planeP - POI, planeN) / dot(lineU, planeN);
    return POI + lineU * d;
}

void main()
{
    // gl_Position = view  * vec4(aPos, 1.0);

    vec4 pos = model * vec4(aPos, 1.0);

    float x = pos.x;
    float z = pos.z;
    float t = uTime;

    Normal = vec3(0, 1, 0);

    float L = 75.0;
    float A = 0.5;
    float S = 40.0;

    float persistance = 0.5;
    float lacunarity = 0.72;
    float prevDx = 0.0;
    float prevDz = 0.0;

    for (float i = 0; i < 6.0; i++){ //random(vec2(x, i))*5, random(vec2(z, i)))*5
        Wave wv = Wave(L, A, S, normalize(
                                    vec2(
                                        random(i),
                                        random(0.5+i)
                                    )
                                )
                            *5.0);
    //     // Wave wv = Wave(L, A, 5, vec2(random(vec2(x, i)), random(vec2(z, i))));
    //     pos.y += waveFunction(x + prevDx, z + prevDz, t, wv);

        prevDx = waveDerivativeX(x, z, t, wv);
        prevDz = waveDerivativeZ(x, z, t, wv);

        Normal.x += prevDx;
        Normal.z += prevDz;


        L *= lacunarity;
        A *= persistance;
        S *= lacunarity;
    }

    // Begin Caustic Calculation
    // Snell's law
    vec3 N = normalize(Normal);
    vec3 ray = normalize(refract(vec3(0.0, 1.0, 0.0), N, 1.0/1.33));
    vec3 POI = vec3(pos.x, 0.0, pos.z);


    oldPos = pos.xyz;
    newPos = line_plane_intercept(POI, pos.xyz, vec3(0.0, 1.0, 0.0), ray);
    // Caustic = abs(asin(sqrt(1.0 - pow(dot(N, vec3(0.0, 1.0, 0.0)), 2.0))/1.33));

    FragPos = pos.xyz/pos.w;
    gl_Position = projection * view * pos;

    Selected = selected;
    TexCoord = vec2(aTexCoord.x, 1.0 - aTexCoord.y);
}
