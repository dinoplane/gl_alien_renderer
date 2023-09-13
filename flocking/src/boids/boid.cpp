#include <glad/glad.h>
#include <GLFW/glfw3.h>

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


std::vector<Boid*>  Boid::boids;
// std::vector<Boid*>  Boid::exper;

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
    // std::cout << "Init Force: " << glm::to_string(glm::vec4(force, 1.0))  << " : " <<  force.length() << std::endl;

    for (auto pair : behaviors){
        force += pair.first * pair.second(*this);
    }
    // force.y = 0.0f;
    // std::cout << "Calc Force: " << glm::to_string(glm::vec4(force, 1.0))  << " : " <<  glm::length(force) << std::endl;

    // if (glm::length(force) > maxForce){
    //     force = glm::normalize(force) * maxForce;
    // }
    // std::cout << "Limi Force: " << glm::to_string(glm::vec4(force, 1.0)) << std::endl;
}

void Boid::updatePosition(){
    acceleration += force;
    // glm::vec3 ne_v = velocity + acceleration;
    // std::cout << "v: " << glm::to_string(velocity) << std::endl;
    // std::cout << "f: " << glm::to_string(force) << std::endl;
    velocity += acceleration;
    // if (glm::normalize(ne_v) != glm::normalize(velocity)){
    //     std::cout << glfwGetTime() << ": Switch: " <<
    //         glm::to_string(velocity) << " -> " <<
    //         glm::to_string(ne_v) << std::endl;
    // }
    // velocity = ne_v;
    // std::cout << "Velocity: " << glm::to_string(glm::vec4(velocity,1.0)) << std::endl;
    // clamp the velocity value by the magnitude of maxSpeed
    velocity = limit(velocity, maxSpeed);

    // glm::vec3 lim = limit(velocity, maxSpeed);
    // if (glm::normalize(lim) != glm::normalize(velocity)){
    //     std::cout << glfwGetTime() << ": BAM: " <<
    //         glm::to_string(velocity) << " -> " <<
    //         glm::to_string(lim) << std::endl;
    // }
    // velocity = lim;

    position += velocity;
    acceleration *= 0.0f;
    force *= 0.0f;
    model.reset();

    model.translate(position);

    model.modelMat *= glm::inverse(glm::lookAt(glm::vec3(0.0), velocity, glm::vec3(0.0, 1.0, 0.0)));

    model.rotate(glm::radians(90.0f), glm::vec3(-1.0, 0.0,0.0));
    model.scale(glm::vec3(0.5f, 0.7f, 0.5f));

    // std::cout << "Position: " << glm::to_string(glm::vec4(position, 1.0)) << std::endl;

    // std::cout << glm::to_string(model.modelMat * glm::vec4(position, 1.0)) << std::endl;
}


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

glm::vec3 Boid::borders(Boid& b){
    glm::vec3 steer = glm::vec3(0.0);
    // if ((abs(b.position.x) > 10) || (abs(b.position.z) > 10) )
    //     std::cout << glfwGetTime() << ": BOING" << std::endl;
    steer.x = (abs(b.position.x) > 10) ? glm::sign(b.position.x) * -0.0005f * (abs(b.position.x) - 10): 0; //
    steer.y = (abs(b.position.y) > 10) ? glm::sign(b.position.y) * -0.0005f * (abs(b.position.y) - 10): 0; //
    steer.z = (abs(b.position.z) > 10) ? glm::sign(b.position.z) * -0.0005f * (abs(b.position.z) - 10): 0; //
    return steer;
}

glm::vec3 Boid::seek(Boid& b, glm::vec3 target){
    glm::vec3 desired = target - b.position;
    desired = glm::normalize(desired);
    desired *= b.maxSpeed;
    desired -= b.velocity;
    glm::vec3 steer = limit(desired, b.maxForce);
    return steer;
}

glm::vec3 Boid::cohesion(Boid& b){
    float neighborDist = 2.5;

    glm::vec3 avgPos = glm::vec3(0.0f);   // Start with empty vector to accumulate all locations
    glm::vec3 steer = glm::vec3(0.0f);
    std::vector<Boid*> flock = Boid::getAllBoidsInRange(b, neighborDist);
    for (Boid* d : flock){
        avgPos += d->position;
    }
    if (flock.size() > 0){
        avgPos /= flock.size();
        steer = Boid::seek(b, avgPos);
    }

    std::cout << "ID" << b.ID <<
    ": Cohesion: " << glm::to_string(steer)
    << " " << flock.size() << " "<<  std::endl;

    return steer;
}


glm::vec3 Boid::separation(Boid& b){
    float desiredSeparation = 1.2;
    glm::vec3 steer = glm::vec3(0.0);
    // std::cout << "copy" << std::endl;
    std::vector<Boid*> flock = Boid::getAllBoidsInRange(b, desiredSeparation);
    glm::vec3 diff;
    // if (flock.size() > 0)
    //     std::cout << b.ID<<  " "<< flock.size() << std::endl;
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

    // std::cout << "ID" << b.ID << ": Separation: " << glm::to_string(steer) << std::endl;

    return steer;
}

glm::vec3 Boid::alignment(Boid& b){
    float neighborDist = 2.0;
    glm::vec3 steer = glm::vec3(0.0);
    std::vector<Boid*> flock = Boid::getAllBoidsInRange(b, neighborDist);
    // if (flock.size() != 0 && flock.size() != 100)
    //     std::cout << flock.size() << std::endl;
    for (Boid* d : flock){
        steer += d->velocity;
    }
    if (flock.size() > 0){
        steer /= flock.size();

        steer = glm::normalize(steer);
        steer *= b.maxSpeed;
        steer -= b.velocity;
    }

    if (glm::length(steer) > b.maxForce){
        steer = glm::normalize(steer) * b.maxForce;
    }
    // std::cout << "Alignment: " << glm::to_string(steer) << std::endl;
    return steer;
}

void Boid::render(){

    this->model.render();

}

std::vector<Boid*> Boid::getAllBoidsInRange(Boid& b, float range){
    //         return Boid::boids | std::ranges::views::filter([](Boid d) {
    //     return glm::length(b.position - d.position) < range;
    // });
    // std::cout << "copy end" << std::endl;

        std::vector<Boid*> ret;
        // Boid::count += 1;

        for(auto d = Boid::boids.begin(); d != Boid::boids.end(); d++){
            float dist = glm::length(b.position - (*d)->position) ;
            // if (b.ID == 0 && d.ID == 0){
            //     std::cout << glm::to_string(b.position) <<
            //     glm::to_string(d.position) << std::endl;
            // }
            // std::cout << "ID" << b.ID << " -> " << (*d)->ID<< ": DIST: " << dist << std::endl;

            // std::cout << glm::to_string(d.position) << std::endl;
            if ((b.ID != (*d)->ID) && 0.0f < dist && dist <= range)
                ret.push_back(*d);
        }
        // if (b.ID == 0 && ret.size() > 0){
        //     std::cout << ret[0].ID << std::endl;
        // }

        return ret;
}


float Boid::calculateKineticEnergy(){
    //         return Boid::boids | std::ranges::views::filter([](Boid d) {
    //     return glm::length(b.position - d.position) < range;
    // });
    float ret = 0.0f;
    for(auto d = Boid::boids.begin(); d != Boid::boids.end(); d++){
        ret += 0.5f * pow(glm::length((*d)->velocity), 2);
    }

    return ret;
}


// std::unique_ptr<Boid> createBoid(glm::vec3 pos  = glm::vec3(0.0f), glm::vec3 vel = glm::vec3(1.0, 0.0, 1.0),){

//     return unique_ptr<>
// }
