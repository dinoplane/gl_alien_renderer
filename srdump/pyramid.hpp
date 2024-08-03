#ifndef PYRAMID_H
#define PYRAMID_H

#include<mesh.hpp>

class Pyramid : public Mesh {
    public:
        virtual void calculateVertices();
};
#endif