#version 460 core
out vec4 FragColor;

layout (location = 10) in vec2 TexCoords;

uniform sampler2D screenTexture;

// const float offset = 1.0 / 300.0;
const vec2 pixels = vec2(180.0, 90.0);

void main(void)
{
  	vec2 p = TexCoords.st;

	p.x -= mod(p.x, 1.0 / pixels.x);
	p.y -= mod(p.y, 1.0 / pixels.y);
    
	vec3 col = texture2D(screenTexture, p).rgb;
	gl_FragColor = vec4(col, 1.0);
}