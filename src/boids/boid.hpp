#ifndef BOID_H
#define BOID_H

#include <glm/glm.hpp>

#include <string>
#include <vector>
#include <unordered_set>

#include <iostream>

#include <memory>

#include <spatial_map.hpp>
#include <pyramid.hpp>


class Boid{
    public:
        static std::vector<Boid*> boids; // consider a concurrent linked list??
        static SpatialMap boid_map;

        static int count;
        float maxForce = 0.001;
        float minForce = 0.000000001;
        float maxSpeed = 0.15;

        int ID;

        glm::vec3 force;
        glm::vec3 acceleration;
        glm::vec3 velocity;
        glm::vec3 position;
        glm::vec3 lastPosition;

        glm::vec3 dimensions = glm::vec3(1.0f);

        std::vector<std::pair<float, glm::vec3(*)(Boid&)>> behaviors; // Consider making this a class to make the mutation easier
        // the genetic multithreading will be a fun multithreading experiment

        Pyramid model;

        // SpatialEntry* entry = nullptr;

    // private:

        // Boid() {model.init();};
        Boid(glm::vec3 pos = glm::vec3(0.0f), glm::vec3 vel = glm::vec3(1.0, 0.0, 1.0)) {
            // std::cout<< "hello" <<std::endl;
            position = pos;
            velocity = vel;
            // std::cout << glm::to_string(position) << std::endl;

            acceleration = glm::vec3(0.0f);
            force = glm::vec3(0.0f);
            model.init();
            behaviors.push_back(std::make_pair(1.0f, borders));
            behaviors.push_back(std::make_pair(3.0f, separation));
            behaviors.push_back(std::make_pair(1.5f, alignment));
            behaviors.push_back(std::make_pair(1.0f, cohesion));

            // std::cout << "DEFAULT CALLED" << std::endl;
            ID = Boid::count;
            addToBoidList();
            Boid::count += 1;
        }

        Boid(const Boid &other){
            std::cout<< "hello" <<std::endl;

            std::cout << "COPY CALLED" << std::endl;
            addToBoidList();
        }


        Boid(Boid &&other){
            std::cout << "MOVE CALLED" << std::endl;
            addToBoidList();
        }

        void addToBoidList(){
            Boid::boids.push_back(this);
            Boid::boid_map.insert(this);
        }

    public:
        static std::unordered_set<Boid*> getAllBoidsInRange(Boid & b, float range);
        // static std::unique_ptr<Boid> createBoid(glm::vec3 pos, glm::vec3 vel);

        static float calculateKineticEnergy();

        // Rule of five O_O
        ~Boid(){
            Boid::boids.erase(
                std::find(Boid::boids.begin(),
                Boid::boids.end(), this));
            Boid::boid_map.remove(this);
        }

        Boid& operator=(const Boid &other){
            std::cout<< "hello" <<std::endl;

            return *this;
        }

        Boid& operator=(Boid &&other){
            std::cout<< "hello" <<std::endl;

            return *this;
        }

        void calculateForce();

        void updatePosition();

        void updateMapEntry();


        void update();

        // Absolutes of life
        // birth()???
        // death()???

        // Behaviors
        static glm::vec3 hunger(Boid& b);
        static glm::vec3 reproduction(Boid& b);
        static glm::vec3 borders(Boid& b);
        static glm::vec3 seek(Boid &b, glm::vec3 target);
        static glm::vec3 cohesion(Boid &b);
        static glm::vec3 separation(Boid &b);
        static glm::vec3 alignment(Boid &b);

        void render();
};
#endif