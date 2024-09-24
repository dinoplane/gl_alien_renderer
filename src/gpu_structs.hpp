#ifndef GPU_STRUCTS_H

#define GPU_STRUCTS_H

struct alignas(16) GPUSphere
{                               // base alignment | aligned offset
    float centerAndRadius[4];   // 4 | 0
                                // 4 | 4
                                // 4 | 8
                                // 4 | 12
}; // 64 bytes???


struct alignas(16) GPUPlane
{                               // base alignment | aligned offset
    float normalAndDistance[4]; // 4 | 0
                                // 4 | 4
                                // 4 | 8
                                // 4 | 12
}; // 64 bytes???

struct alignas(16) GPUFrustum
{
    GPUPlane topFace;
    GPUPlane bottomFace;

    GPUPlane rightFace;
    GPUPlane leftFace;

    GPUPlane farFace;
    GPUPlane nearFace;
};

#endif
