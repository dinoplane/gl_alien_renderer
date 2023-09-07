#ifndef SOLID_H
#define SOLID_H
#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

class Solid {
    protected:
        // Create a renderer class and compare times

        unsigned int VBO, VAO, EBO;
        std::vector<float> vertices;
        std::vector<unsigned int> indices;
        // std::vector<float> tricolors;
        virtual void calculateVertices();
        void initializeBuffers();

        virtual void init();

    public:
        std::string type;
        glm::vec4 color;
        glm::mat4 model;

    // private:
    //     void drawTriangle3D();
    //     void pushTriangle3D(float verts[]);


        Solid();
        void deleteBuffers();



        void translate(glm::vec3 delta){
            model = glm::translate(model, delta);
        };
        void rotate(float angle, glm::vec3 axis){
            model = glm::rotate(model, angle, axis);
        };
        void scale(glm::vec3 factor){
            model = glm::scale(model, factor);
        };

        void reset() {


        };

        void render(); // I believe we can optimize this from asg3
};
#endif