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
        glm::mat4 modelMat;

    // private:
    //     void drawTriangle3D();
    //     void pushTriangle3D(float verts[]);


        Solid();
        void deleteBuffers();



        void translate(glm::vec3 delta){
            modelMat = glm::translate(modelMat, delta);
        };
        void rotate(float angle, glm::vec3 axis){
            modelMat = glm::rotate(modelMat, angle, axis);
        };
        void scale(glm::vec3 factor){
            modelMat = glm::scale(modelMat, factor);
        };

        void reset() {
            modelMat = glm::mat4(1.0);

        };

        virtual void render(); // I believe we can optimize this from asg3
};
#endif