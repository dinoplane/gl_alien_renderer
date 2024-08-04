#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


struct Transform {
    glm::vec3 position;
    glm::vec3 rotation; // euler angles
    glm::vec3 scale;

    bool cached = false;
    glm::mat4 modelMat;

    glm::mat4 GetModelMatrix(){
        if(!cached){
            glm::mat4 tmpModelMat = glm::identity<glm::mat4>();
            tmpModelMat = glm::scale(tmpModelMat, scale);
            tmpModelMat = glm::rotate(tmpModelMat, glm::radians(rotation.x), glm::vec3(1.0, 0.0, 0.0));
            tmpModelMat = glm::rotate(tmpModelMat, glm::radians(rotation.y), glm::vec3(0.0, 1.0, 0.0));
            tmpModelMat = glm::rotate(tmpModelMat, glm::radians(rotation.z), glm::vec3(0.0, 0.0, 1.0));
            tmpModelMat = glm::translate(tmpModelMat, position);
            modelMat = tmpModelMat;
            cached = true;
        }
        return modelMat;
    }

    void SetPosition(const glm::vec3& pos){
        position = pos;
        cached = false;
    }

    void SetRotation(const glm::vec3& rot){
        rotation = rot;
        cached = false;
    }

    void SetScale(const glm::vec3& scl){
        scale = scl;
        cached = false;
    }
};

#endif