float sdCircle( vec2 p, float r )
{
    return length(p) - r;
}

void mainImage( out vec4 fragColor, in vec2 fragCoord )
{

    // Normalized pixel coordinates (from 0 to 1)
    vec2 uv = fragCoord/iResolution.xy * 2.0 - 1.0;
    uv.x *=  iResolution.x/iResolution.y;

    // Time varying pixel color
    vec3 col = vec3(0.0, 0.9, 0.9);

    // float d = 1./abs(sdCircle(uv, 0.2));
    float d = sdCircle(uv, 0.5);

    col = mix(vec3(0.0), col, smoothstep(-0.2, 0.1, -abs(d)));
    if (d < 0.0)
        col = mix(vec3(1.0), col, smoothstep(2., -1., d));
    // col = mix(col, vec3(0.0), 1.0-smoothstep(0.0,0.01,abs(d)) );
    // smoothstep(2, -1, x)
    // Output to screen
    // if (d < )
    fragColor = vec4(col,1.0);
}