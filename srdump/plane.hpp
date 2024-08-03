#ifndef PLANE_H
#define PLANE_H

#include <glm/glm.hpp>

#include<mesh.hpp>


class Plane : public Mesh {
    public:

        glm::ivec2 grid_dims = glm::ivec2(1);
        glm::vec2 unit_dims = glm::vec2(1.0f);

        Plane() = delete;
        Plane(glm::ivec2 gdims);

        // virtual void init();
        virtual void calculateVertices();

        void setTiles(glm::ivec2 gdims);
};
#endif