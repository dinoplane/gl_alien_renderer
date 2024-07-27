#ifndef PYRAMID_H
#define PYRAMID_H

#include<solid.hpp>

class Pyramid : public Solid {
    public:
        virtual void init();
        virtual void calculateVertices();
};
#endif