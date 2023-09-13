#ifndef CAMERA_H
#define CAMERA_H

#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


#include <string>
#include <vector>

#include <fstream>
#include <sstream>
#include <iostream>

const float YAW         = -90.0f;
const float PITCH       =  -89.0f;
const float SPEED       =  10.0f;
const float SENSITIVITY =  0.1f;
const float ZOOM        =  45.0f;


class Camera {
    // float fov;
    float speed;
    float sensitivity;

    glm::vec3 front;
    glm::vec3 right;
    glm::vec3 position;
    glm::vec3 worldUp;
    glm::vec3 up;

    // glm::mat4 viewMat;

    glm::mat4 projMat;


    float yaw;
    float pitch;


    public:
        Camera(
                float width, float height,
                glm::vec3 pos = glm::vec3(0.0, 20.0f, 0.0),
                glm::vec3 upv = glm::vec3(0.0, 1.0, 0.0),
                float yaw_val = YAW, float pitch_val = PITCH
            ) : front(glm::vec3(0.0f, 0.0f, -1.0f)), speed(SPEED), sensitivity(SENSITIVITY){
            // viewMat = glm::mat4(1.0);
            position = pos;
            worldUp = upv;
            yaw = yaw_val;
            pitch = pitch_val;
            updateCameraVectors();

            projMat = glm::perspective(glm::radians(60.0f), width / height, 0.1f, 2000.0f);
        }

        void moveLeft(float deltaTime);
        void moveRight(float deltaTime);
        void moveForward(float deltaTime);
        void moveBackward(float deltaTime);

        void processMouseMovement(float xoffset, float yoffset);


        glm::mat4 getViewMatrix(){
            return glm::lookAt(position, position + front, up);
        }

        glm::mat4 getProjMatrix(){
            return projMat;
        }

    private:
        void updateCameraVectors();

};
#endif