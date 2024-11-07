#ifndef COMP_RENDERER_H
#define COMP_RENDERER_H


#include <renderer.hpp>
#include <shader_c.hpp>


class ComputeRenderer : public Renderer
{
    public:
    const unsigned int TEXTURE_WIDTH = 512, TEXTURE_HEIGHT = 512;
    ComputeShader* computeShader;
    uint32_t outTexture;

    ComputeRenderer(float w, float h);

    void Render(float currentFrame);
};

#endif