#version 460 core
out vec4 FragColor;

in vec2 TexCoords;

uniform sampler2D screenTexture;

const float time = 0.5;
const vec2 mouse = vec2(0.0, 0.0);
const float turns = 2.5;

void main(void) {
  vec2 p = -1.0 + 2.0 * TexCoords.st;
  vec2 m = -1.0 + 2.0 * mouse.xy;

  float a1 = abs(atan(p.y - m.y, p.x - m.x));
  float r1 = sqrt(dot(p - m, p - m));
  float a2 = abs(atan(p.y + m.y, p.x + m.x));
  float r2 = sqrt(dot(p + m, p + m));

  vec2 uv;
  uv.x = time + (r1 - r2) * 0.25;
  uv.y = sin(turns * (a1 - a2));

  float w = r1 * r2 * 0.5;
  vec3 col = texture2D(screenTexture, 0.5 - 0.495 * uv).xyz;
  gl_FragColor = vec4(col / (0.01 + w), 1.0);
}