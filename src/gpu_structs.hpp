#ifndef GPU_STRUCTS_H

#define GPU_STRUCTS_H

struct GPUSphere
{
    float centerAndRadius[4];
};


struct GPUPlane
{
    float normalAndDistance[4];
};

struct GPUFrustum
{
    GPUPlane topFace;
    GPUPlane bottomFace;

    GPUPlane rightFace;
    GPUPlane leftFace;

    GPUPlane farFace;
    GPUPlane nearFace;
};

#endif
