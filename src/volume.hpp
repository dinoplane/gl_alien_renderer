#ifndef VOLUME_H
#define VOLUME_H


#include <util.h>
#include <transform.hpp>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


#include <string>
#include <iostream>
#include <gpu_structs.hpp>


struct Volume
{

    virtual bool IsOnFrustum(const Frustum& camFrustum,
                                            const Transform& modelTransform) const = 0;
};

struct Sphere : public Volume
{
    glm::vec3 center{ 0.f, 0.f, 0.f };
    float radius{ 0.f };

    Sphere()
        : center{ 0.f, 0.f, 0.f }, radius{ 0.f }
    {}

    Sphere(const glm::vec3& inCenter, float inRadius)
        : center{ inCenter }, radius{ inRadius }
    {}

    bool IsOnOrForwardPlane(const Plane& plane) const
    {
        return plane.GetSignedDistanceToPlane(center) > -radius;
    }

    bool IsOnFrustum(const Frustum& camFrustum,
                                            const Transform& modelTransform) const final;

    GPUSphere ToGPUSphere() const 
    {
        return { {center.x, center.y, center.z, radius} };
    }
};
/*
struct AABB : public BoundingVolume
{
    glm::vec3 center{ 0.f, 0.f, 0.f };
    glm::vec3 extents{ 0.f, 0.f, 0.f };

    AABB(const glm::vec3& min, const glm::vec3& max)
        : BoundingVolume{},
        center{ (max + min) * 0.5f },
        extents{ max.x - center.x, max.y - center.y, max.z - center.z }
    {}

    AABB(const glm::vec3& inCenter, float iI, float iJ, float iK)
        : BoundingVolume{}, center{ inCenter }, extents{ iI, iJ, iK }
    {}

    bool isOnOrForwardPlane(const Plane& plane) const
    {
        // Compute the projection interval radius of b onto L(t) = b.c + t * p.n
        const float r = extents.x * std::abs(plane.normal.x) +
                extents.y * std::abs(plane.normal.y) + extents.z * std::abs(plane.normal.z);

        return -r <= plane.getSignedDistanceToPlane(center);
    }

    bool isOnFrustum(const Frustum& camFrustum, const Transform& transform) const final
    {
        //Get global scale thanks to our transform
        const glm::vec3 globalCenter{ transform.GetModelMatrix() * glm::vec4(center, 1.f) };

        // Scaled orientation
        const glm::vec3 right = transform.getRight() * extents.x;
        const glm::vec3 up = transform.getUp() * extents.y;
        const glm::vec3 forward = transform.getForward() * extents.z;

        const float newIi = std::abs(glm::dot(glm::vec3{ 1.f, 0.f, 0.f }, right)) +
            std::abs(glm::dot(glm::vec3{ 1.f, 0.f, 0.f }, up)) +
            std::abs(glm::dot(glm::vec3{ 1.f, 0.f, 0.f }, forward));

        const float newIj = std::abs(glm::dot(glm::vec3{ 0.f, 1.f, 0.f }, right)) +
            std::abs(glm::dot(glm::vec3{ 0.f, 1.f, 0.f }, up)) +
            std::abs(glm::dot(glm::vec3{ 0.f, 1.f, 0.f }, forward));

        const float newIk = std::abs(glm::dot(glm::vec3{ 0.f, 0.f, 1.f }, right)) +
            std::abs(glm::dot(glm::vec3{ 0.f, 0.f, 1.f }, up)) +
            std::abs(glm::dot(glm::vec3{ 0.f, 0.f, 1.f }, forward));

        //We not need to divise scale because it's based on the half extention of the AABB
        const AABB globalAABB(globalCenter, newIi, newIj, newIk);

        return (globalAABB.isOnOrForwardPlane(camFrustum.leftFace) &&
            globalAABB.isOnOrForwardPlane(camFrustum.rightFace) &&
            globalAABB.isOnOrForwardPlane(camFrustum.topFace) &&
            globalAABB.isOnOrForwardPlane(camFrustum.bottomFace) &&
            globalAABB.isOnOrForwardPlane(camFrustum.nearFace) &&
            globalAABB.isOnOrForwardPlane(camFrustum.farFace));
    };
};

*/
#endif