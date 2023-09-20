#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#include<boid.hpp>
#include <pyramid.hpp>
// #include <spatial_map.hpp>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>


std::vector<Boid*>  Boid::boids;
int Boid::count = 0;

glm::vec3 limit(glm::vec3 v, float maxMag){
    glm::vec3 ret = v;
    if (glm::length(v) > maxMag){
        ret = glm::normalize(v) * maxMag;
    }
    return ret;
}

void Boid::calculateForce(){
    force = glm::vec3(0.0);

    for (auto pair : behaviors){
        force += pair.first * pair.second(*this);
    }

    // std::string msg = std::to_string(ID) + "  finished\n";
    // std::cout << msg;
}

void Boid::updatePosition(){
    // Update acceleration and velocity with cached force
    acceleration += force;
    velocity += acceleration;
    velocity = limit(velocity, maxSpeed);

    // Update position with velocity
    position += velocity;

    // Reset acceleration and force
    acceleration *= 0.0f;
    force *= 0.0f;

    // Apply necessary transforms to model
    model.reset();
    model.translate(position);
    model.modelMat *= glm::inverse(glm::lookAt(glm::vec3(0.0), velocity, glm::vec3(0.0, 1.0, 0.0)));
    model.rotate(glm::radians(90.0f), glm::vec3(-1.0, 0.0,0.0));
    model.scale(glm::vec3(0.5f, 0.7f, 0.5f));

    // std::cout << "Position: " << glm::to_string(glm::vec4(position, 1.0)) << std::endl;
}

// Not called cuz multithreading
void Boid::update(){
    calculateForce();
    updatePosition();
}

glm::vec3 Boid::hunger(Boid& b){
    return glm::vec3(0.0f);
}

glm::vec3 Boid::reproduction(Boid& b){
    return glm::vec3(0.0f);
}

/*
    Apply force to boids proportional to distance outside the boundary
*/
glm::vec3 Boid::borders(Boid& b){
    glm::vec3 steer = glm::vec3(0.0);

    steer.x = (abs(b.position.x) > 10) ? glm::sign(b.position.x) * -0.0005f * (abs(b.position.x) - 10): 0; //
    steer.y = (abs(b.position.y) > 10) ? glm::sign(b.position.y) * -0.0005f * (abs(b.position.y) - 10): 0; //
    steer.z = (abs(b.position.z) > 10) ? glm::sign(b.position.z) * -0.0005f * (abs(b.position.z) - 10): 0; //

    return steer;
}


/*
    Seek towards target.
*/
glm::vec3 Boid::seek(Boid& b, glm::vec3 target){
    glm::vec3 desired = target - b.position;
    desired = glm::normalize(desired);
    desired *= b.maxSpeed;
    desired -= b.velocity;
    glm::vec3 steer = limit(desired, b.maxForce);

    return steer;
}

/*
    Mingle with your neighbors.
*/
glm::vec3 Boid::cohesion(Boid& b){
    float neighborDist = 2.5;
    glm::vec3 avgPos;   // Start with empty vector to accumulate all locations
    glm::vec3 steer = glm::vec3(0.0f);
    std::vector<Boid*> flock = Boid::getAllBoidsInRange(b, neighborDist);

    for (Boid* d : flock){
        avgPos += d->position;
    }

    if (flock.size() > 0){
        avgPos /= flock.size();
        steer = Boid::seek(b, avgPos);
    }

    return steer;
}

/*
    Move away from your closest neighbors
*/
glm::vec3 Boid::separation(Boid& b){
    float desiredSeparation = 1.2;
    glm::vec3 steer = glm::vec3(0.0);
    std::vector<Boid*> flock = Boid::getAllBoidsInRange(b, desiredSeparation);
    glm::vec3 diff;

    for (Boid* d : flock){
        diff = b.position - d->position;
        float dist = glm::length(diff);
        diff = glm::normalize(diff);
        diff /= (dist*dist);
        steer += diff;
    }

    if (flock.size() > 0)
        steer /= flock.size();

    if (glm::length(steer) > 0.0000001){
        steer = glm::normalize(steer);
        steer *= b.maxSpeed;
        steer -= b.velocity;
        steer = limit(steer, b.maxForce);
    }

    return steer;
}

/*
    Match your direction with your neighbors
*/
glm::vec3 Boid::alignment(Boid& b){
    float neighborDist = 2.0;
    glm::vec3 steer = glm::vec3(0.0);
    std::vector<Boid*> flock = Boid::getAllBoidsInRange(b, neighborDist);

    for (Boid* d : flock){
        steer += d->velocity;
    }

    if (flock.size() > 0){
        steer /= flock.size();

        steer = glm::normalize(steer);
        steer *= b.maxSpeed;
        steer -= b.velocity;
        steer = limit(steer, b.maxForce);
    }

    return steer;
}

void Boid::render(){

    this->model.render();

}

std::vector<Boid*> Boid::getAllBoidsInRange(Boid& b, float range){
    std::vector<Boid*> ret;

    for(auto d = Boid::boids.begin(); d != Boid::boids.end(); d++){
        float dist = glm::length(b.position - (*d)->position) ;

        if ((b.ID != (*d)->ID) && 0.0f < dist && dist <= range)
            ret.push_back(*d);
    }

    return ret;
}


float Boid::calculateKineticEnergy(){
    float ret = 0.0f;
    for(auto d = Boid::boids.begin(); d != Boid::boids.end(); d++){
        ret += 0.5f * pow(glm::length((*d)->velocity), 2);
    }

    return ret;
}


// std::unique_ptr<Boid> createBoid(glm::vec3 pos  = glm::vec3(0.0f), glm::vec3 vel = glm::vec3(1.0, 0.0, 1.0),){

//     return unique_ptr<>
// }
