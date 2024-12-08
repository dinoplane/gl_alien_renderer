#version 460 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture;
const float abberationAmount = 0.01; // Control the intensity of the aberration

void main() {
    // Offset for chromatic aberration
    vec2 redOffset = vec2(abberationAmount*sqrt(3)*0.5, -0.5*abberationAmount);
    vec2 greenOffset = vec2(-abberationAmount*sqrt(3)*0.5, -0.5*abberationAmount);
    vec2 blueOffset = vec2(0.0, abberationAmount);

    // Sample the texture at different offsets for RGB channels
    float r = texture(screenTexture, TexCoords + redOffset).r;
    float g = texture(screenTexture, TexCoords + greenOffset).g;
    float b = texture(screenTexture, TexCoords + blueOffset).b;

    // Combine channels into final color
    FragColor = vec4(r, g, b, 1.0);
}