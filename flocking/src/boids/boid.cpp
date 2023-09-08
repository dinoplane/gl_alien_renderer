#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <glm/gtx/string_cast.hpp>

#include<boid.hpp>
#include <pyramid.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

void Boid::calculateForce(){
    force = glm::vec3(0.0);
    for (auto pair : behaviors){
        force += pair.first * pair.second();
    }
}

void Boid::updatePosition(){
    this->acceleration += this->force;
    this->velocity += this->acceleration;
    // clamp the velocity value by the magnitude of maxSpeed
    if (this->velocity.length() > this->maxSpeed){
        this->velocity = glm::normalize(this->velocity) * this->maxSpeed;
    }

    this->position += this->velocity;
    this->acceleration *= 0.0f;
    this->model.reset();

    this->model.translate(this->position);

    this->model.modelMat *= glm::lookAt(glm::vec3(0.0), this->velocity, glm::vec3(0.0, 1.0, 0.0));

    this->model.rotate(glm::radians(90.0f), glm::vec3(0.0, 0.0,1.0));
    this->model.scale(glm::vec3(0.5f, 0.7f, 0.5f));

    // std::cout << "Position: " << glm::to_string(glm::vec4(position, 1.0)) << std::endl;

    // std::cout << glm::to_string(model.modelMat * glm::vec4(position, 1.0)) << std::endl;
}

glm::vec3 Boid::hunger(){
    return glm::vec3(0.0f);
}

glm::vec3 Boid::reproduction(){
    return glm::vec3(0.0f);
}

glm::vec3 Boid::borders(){
    glm::vec3 steer = glm::vec3(0.0);
    steer.x = (abs(position.x) > 100) ? glm::sign(position.x) * -25 : 0;
    steer.y = (abs(position.y) > 100) ? glm::sign(position.y) * -25 : 0;
    steer.z = (abs(position.z) > 100) ? glm::sign(position.z) * -25 : 0;
    return steer;
}

glm::vec3 Boid::cohesion(){
    return glm::vec3(0.0f);
}
glm::vec3 Boid::avoidance(){
    return glm::vec3(0.0f);
}

glm::vec3 Boid::alignment(){
    return glm::vec3(0.0f);
}

void Boid::render(){

    this->model.render();

}

