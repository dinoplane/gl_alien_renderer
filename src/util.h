#ifndef UTIL_H
#define UTIL_H

#include <glm/glm.hpp>
#include <gpu_structs.hpp>

#define ARRAY_COUNT(x) (sizeof(x) / sizeof((x)[0]))

template <typename T, size_t N>
constexpr size_t array_size(T (&)[N]) {
    return N;
}


struct Vertex{
    glm::vec3 position;
    glm::vec3 normal;
    glm::vec2 texcoords;
};




struct Plane
{
    // unit vector
    glm::vec3 normal = { 0.f, 1.f, 0.f };

    // distance from origin to the nearest point in the plane
    float     distance = 0.f;

	Plane() = default;

	Plane(const glm::vec3& p1, const glm::vec3& norm)
		: normal(glm::normalize(norm)),
		distance(glm::dot(normal, p1))
        // the normal is parallel to the vector from the nearest point to the origin
        // u dot v is measures the component of u in the v direction (where the projection equation comes from)
	{}

	float GetSignedDistanceToPlane(const glm::vec3& point) const
	{
		return glm::dot(normal, point) - distance;
	}

    GPUPlane ToGPUPlane()
    {
        return { {normal.x, normal.y, normal.z, distance} };
    }
};

struct Frustum
{
    Plane topFace;
    Plane bottomFace;

    Plane rightFace;
    Plane leftFace;

    Plane farFace;
    Plane nearFace;

    GPUFrustum ToGPUFrustum()
    {
        return { topFace.ToGPUPlane(), bottomFace.ToGPUPlane(),
            rightFace.ToGPUPlane(), leftFace.ToGPUPlane(),
            farFace.ToGPUPlane(), nearFace.ToGPUPlane() };
    }
};


#endif