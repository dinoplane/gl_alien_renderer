#ifndef ENTITY_H
#define ENTITY_H

#include "model.hpp"
#include "transform.hpp"


struct Entity {
    Model model;
    Transform transform;
    // Add a material later

    Entity(const Model& inModel, const Transform& inTransform)
        : model{ inModel }, transform{ inTransform }
    {}
};

#endif