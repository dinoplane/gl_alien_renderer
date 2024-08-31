
#include <volume.hpp>
#include <util.h>

#include <transform.hpp>
#include <iostream>


bool Sphere::IsOnFrustum(const Frustum& camFrustum,
                                            const Transform& modelTransform) const
{
    //Get global scale is computed by doing the magnitude of
    //X, Y and Z model matrix's column.
    const glm::vec3 globalScale = modelTransform.GetGlobalScale();

    //Get our global center with process it with the global model matrix of our modelTransform
    const glm::vec3 globalCenter (modelTransform.GetModelMatrix() * glm::vec4(center, 1.f));

    //To wrap correctly our shape, we need the maximum scale scalar.
    const float maxScale = std::max(std::max(globalScale.x, globalScale.y), globalScale.z);

    Sphere globalSphere(globalCenter, radius * maxScale);

    //Check Firstly the result that have the most chance
    //to faillure to avoid to call all functions.
    return (globalSphere.IsOnOrForwardPlane(camFrustum.leftFace) &&
        globalSphere.IsOnOrForwardPlane(camFrustum.rightFace) &&
        globalSphere.IsOnOrForwardPlane(camFrustum.farFace) &&
        globalSphere.IsOnOrForwardPlane(camFrustum.nearFace) &&
        globalSphere.IsOnOrForwardPlane(camFrustum.topFace) &&
        globalSphere.IsOnOrForwardPlane(camFrustum.bottomFace));
};