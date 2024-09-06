#ifndef TRANSFORM_H
#define TRANSFORM_H

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>


struct Transform {
    glm::vec3 position = glm::vec3(0.0, 0.0, 0.0);
    glm::vec3 rotation = glm::vec3(0.0, 0.0, 0.0); // euler angles
    glm::vec3 scale = glm::vec3(1.0, 1.0, 1.0);

    glm::mat4 GetModelMatrix() const {

        glm::mat4 tmpModelMat = glm::identity<glm::mat4>();
        tmpModelMat = glm::scale(tmpModelMat, scale);
        tmpModelMat = glm::rotate(tmpModelMat, glm::radians(rotation.x), glm::vec3(1.0, 0.0, 0.0));
        tmpModelMat = glm::rotate(tmpModelMat, glm::radians(rotation.y), glm::vec3(0.0, 1.0, 0.0));
        tmpModelMat = glm::rotate(tmpModelMat, glm::radians(rotation.z), glm::vec3(0.0, 0.0, 1.0));
        tmpModelMat = glm::translate(tmpModelMat, position);

        return tmpModelMat;
    }

    glm::vec3 GetGlobalScale() const {
        return scale;
    }

    void SetPosition(const glm::vec3& pos){
        position = pos;
    }

    void SetRotation(const glm::vec3& rot){
        rotation = rot;
    }

    void SetScale(const glm::vec3& scl){
        scale = scl;
    }
};

#endif