#include <glad/glad.h>
#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include<camera.hpp>

#include <string>
#include <vector>

#include <fstream>
#include <sstream>
#include <iostream>


// const float SPEED_MULTIPLIER = 10;

// Camera::Camera(float width, float height,
//                 glm::vec3 pos = glm::vec3(0.0, 0.0, 0.0), glm::vec3 upv = glm::vec3(0.0, 1.0, 0.0),
//                 float yaw_val = YAW, float pitch_val = PITCH) :
//                 front(glm::vec3(0.0f, 0.0f, -1.0f)), speed(SPEED), sensitivity(SENSITIVITY){
//     // viewMat = glm::mat4(1.0);
//     position = pos;
//     worldUp = upv;
//     yaw = yaw_val;
//     pitch = pitch_val;
//     updateCameraVectors();

//     projMat = glm::perspective(glm::radians(60.0f), width / height, 0.1f, 2000.0f);
// }

void Camera::moveLeft(float deltaTime){
    position -= right * deltaTime * SPEED;
}

void Camera::moveRight(float deltaTime){
    position += right * deltaTime * SPEED;
}

void Camera::moveForward(float deltaTime){
    position += front * deltaTime * SPEED;
}

void Camera::moveBackward(float deltaTime){
    position -= front * deltaTime * SPEED;
}

void Camera::processMouseMovement(float xoffset, float yoffset){
    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw   += xoffset;
    pitch += yoffset;

    if(pitch > 89.0f)
        pitch = 89.0f;
    if(pitch < -89.0f)
        pitch = -89.0f;

    updateCameraVectors();
}


void Camera::updateCameraVectors(){
    glm::vec3 direction;
    direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    direction.y = sin(glm::radians(pitch));
    direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    front = glm::normalize(direction);
    right = glm::normalize(glm::cross(front, worldUp));
    up = glm::normalize(glm::cross(right, front));
}