#ifndef CUBE_H
#define CUBE_H

#include<mesh.hpp>

class Cube : public Mesh {
    public:
        virtual void calculateVertices();
};
#endif