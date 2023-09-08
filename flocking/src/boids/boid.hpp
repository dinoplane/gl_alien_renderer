#ifndef BOID_H
#define BOID_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include<pyramid.hpp>

#include <string>
#include <vector>

#include <fstream>
#include <sstream>
#include <iostream>

class Boid {

    public:
        static std::vector<Boid> boids; // consider a concurrent linked list??
        float maxForce = 0.05;
        float maxSpeed = 0.01;

        glm::vec3 force;
        glm::vec3 acceleration;
        glm::vec3 velocity;
        glm::vec3 position;
        glm::vec3 cachedPosition;


        std::vector<std::pair<float, glm::vec3(*)()>> behaviors; // Consider making this a class to make the mutation easier
        // the genetic multithreading will be a fun multithreading experiment

        Pyramid model;
        // glm::mat4 modelMat;

    public:
        std::vector<Boid> getAllBoidsInRange(float range);

        // Boid() {model.init();};
        Boid(glm::vec3 pos = glm::vec3(0.0f), glm::vec3 vel = glm::vec3(1.0, -1.0, 1.0)) {
            position = pos;
            velocity = vel;
            model.init();
            updatePosition();
        }

        void calculateForce();

        void updatePosition();

        // Absolutes of life
        // birth()???
        // death()???

        // Behaviors
        glm::vec3 hunger();
        glm::vec3 reproduction();
        glm::vec3 borders();
        glm::vec3 cohesion();
        glm::vec3 avoidance();
        glm::vec3 alignment();

        void render();
};
#endif