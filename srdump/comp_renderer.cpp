#include <comp_renderer.hpp>

ComputeRenderer::ComputeRenderer(float w=512, float h=512) : Renderer(w, h)
{
    computeShader = new ComputeShader("./resources/shader/comp_shader.comp");


    glCreateTextures(GL_TEXTURE_2D, 1, &outTexture);
    glTextureParameteri(outTexture, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTextureParameteri(outTexture, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTextureParameteri(outTexture, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTextureParameteri(outTexture, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTextureStorage2D(outTexture, 1, GL_RGBA32F, TEXTURE_WIDTH, TEXTURE_HEIGHT);
}


void ComputeRenderer::Render(float currentFrame)
{
    glBindImageTexture(0, outTexture, 0, GL_FALSE, 0, GL_READ_WRITE, GL_RGBA32F);
    computeShader->use();
    computeShader->setFloat("t", currentFrame);
    glDispatchCompute((unsigned int) TEXTURE_WIDTH, (unsigned int) TEXTURE_HEIGHT, 1);

    // make sure writing to image has finished before read
    glMemoryBarrier(GL_SHADER_IMAGE_ACCESS_BARRIER_BIT);

    passthroughShader->use();
    glBindVertexArray(Renderer::quadVAO);
    glVertexArrayVertexBuffer(Renderer::quadVAO, 0, Renderer::quadVBO, 0, sizeof(float) * 4);
    glDisable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, outTexture);

    glBindFramebuffer(GL_FRAMEBUFFER, FBO);
    glNamedFramebufferTexture(FBO, GL_COLOR_ATTACHMENT0, dstFBOTexture, 0);
    glDrawArrays(GL_TRIANGLES, 0, 6);
}