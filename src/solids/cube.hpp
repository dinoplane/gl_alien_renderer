#ifndef CUBE_H
#define CUBE_H

#include<solid.hpp>

class Cube : public Solid {
    public:
        virtual void init();
        virtual void calculateVertices();
};
#endif