#ifndef ENTITY_H
#define ENTITY_H

#include "mesh.hpp"
#include "transform.hpp"


struct Entity {
    Mesh mesh;
    Transform transform;
    // Add a material later
};

#endif