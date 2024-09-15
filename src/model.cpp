#include <model.hpp>


Model Model::CreateCube(){
    Model retModel;
    retModel.meshes.push_back(Mesh::CreateCube());
    retModel.nodes.push_back( {glm::identity<glm::mat4x4>(), 0} );
    retModel.boundingVolume = retModel.meshes[0].boundingVolume;
    return retModel;
}

Model Model::CreatePyramid(){
    Model retModel;
    retModel.meshes.push_back(Mesh::CreatePyramid());
    retModel.nodes.push_back( {glm::identity<glm::mat4x4>(), 0} );
    retModel.boundingVolume = retModel.meshes[0].boundingVolume;
    return retModel;
}
