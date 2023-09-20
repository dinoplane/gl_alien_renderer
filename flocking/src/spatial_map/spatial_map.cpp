#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>

#include <unordered_map>
#include <vector>
#include <set> // consider and benchmark unordered
#include <memory>

#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <cstdio>


#include <boid.hpp>
#include <spatial_map.hpp>

/*
    Returns the index of the cell the position is in.

*/
glm::vec3 SpatialMap::_key(glm::vec3 pos){
    std::cout << "Pos - Minbound: "<<  glm::to_string(pos - minbound) << std::endl;
    std::cout << "Bound Dims: "<<  glm::to_string(bound_dims) << std::endl;
    std::cout << "Grid Dims: "<<  glm::to_string(grid_dims) << std::endl;

    return  glm::floor((pos - minbound) / bound_dims * (grid_dims));
}

void SpatialMap::printMap(){
    for (int j = 0; j < grid_dims.y; j++){
        std::cout << "Y = " << j << std::endl;
        std::cout << "   ";
        for (int i = 0; i < grid_dims.x; i++){
            std::printf("  %d  ", i);
        }
        std::cout << std::endl;
        for (int k = 0; k < grid_dims.z; k++){
            std::cout << k << "  ";
            for (int i = 0; i < grid_dims.x; i++){
                glm::vec3 key = glm::vec3(i, j, k);
                std::printf("[ %lu ]", sp_map[key].size());
                // std::unordered_map<glm::vec3, int> sp_set;
                // sp_set[key] = 1;

                // std::cout << glm::to_string(key) << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }

    std::cout << "Out of Bounds: " << sp_map[glm::vec3(-1.0f)].size() << std::endl;
}


void SpatialMap::insert(std::shared_ptr<Boid> b){
    //Lazy insertion? or out of bounds bucket?

    std::cout << "Hello" << std::endl;
    glm::vec3 blb = b->position - b->dimensions/2.0f;
    glm::vec3 trt = b->position + b->dimensions/2.0f;

    glm::vec3 blb_ind = _key(blb);
    glm::vec3 trt_ind = _key(trt);

    std::cout  << "Bottom Left Back: " << glm::to_string(blb_ind) << std::endl;
    std::cout << "Top Right Top: " << glm::to_string(trt_ind) << std::endl;

    // Everything below can be private?
    std::shared_ptr<SpatialEntry> entry(new SpatialEntry(b, b->position, b->dimensions, blb_ind, trt_ind));

    _insert(entry);
}

void SpatialMap::_insert(std::shared_ptr<SpatialEntry> e){
    glm::vec3 blb_ind = e->blb_ind;
    glm::vec3 trt_ind = e->trt_int;

    // bool blb_out = glm::any(glm::lessThanEqual(blb_ind, glm::vec3(0.0f)));
    // bool trt_out =glm::any(glm::greaterThanEqual(trt_ind, grid_dims));

    // if (blb_out && trt_out){
    //     sp_map[glm::vec3(-1.0f)].emplace(b);
    //     blb_ind = glm::vec3(0.0f);
    //     trt_ind = (glm::length(trt - blb) > glm::length(unit_dims)) ? grid_dims - glm::vec3(1.0f); : glm::vec3(0.0f);
    // } else if (blb_out){
    //     sp_map[glm::vec3(-1.0f)].emplace(b);
    //     blb_ind = glm::vec3(0.0f);
    // } else if (trt_out){
    //     sp_map[glm::vec3(-1.0f)].emplace(b);
    //     trt_ind = grid_dims - glm::vec3(1.0f);
    // }


    for (float i = blb_ind.x; i <= trt_ind.x; i++){
        for (float j = blb_ind.y; j <= trt_ind.y; j++){
            for (float k = blb_ind.z; k <= trt_ind.z; k++){
                sp_map[glm::vec3(i, j, k)].emplace(e);
            }
        }
    }


}


void SpatialMap::remove(std::shared_ptr<SpatialEntry> e){
    glm::vec3 blb_ind = e->blb_ind;
    glm::vec3 trt_ind = e->trt_int;

    std::cout  << "Bottom Left Back: " << glm::to_string(blb_ind) << std::endl;
    std::cout << "Top Right Top: " << glm::to_string(trt_ind) << std::endl;

    // bool blb_out = glm::any(glm::lessThanEqual(blb_ind, glm::vec3(0.0f)));
    // bool trt_out =glm::any(glm::greaterThanEqual(trt_ind, grid_dims));

    // if (blb_out && trt_out){
    //     sp_map[glm::vec3(-1.0f)].erase(b);
    //     blb_ind = glm::vec3(0.0f);
    //     trt_ind = (glm::length(trt - blb) > glm::length(unit_dims)) ? grid_dims - glm::vec3(1.0f); : glm::vec3(0.0f);
    // } else if (blb_out){
    //     sp_map[glm::vec3(-1.0f)].erase(b);
    //     blb_ind = glm::vec3(0.0f);
    // } else if (trt_out){
    //     sp_map[glm::vec3(-1.0f)].erase(b);
    //     trt_ind = grid_dims - glm::vec3(1.0f);
    // }


    for (float i = blb_ind.x; i <= trt_ind.x; i++){
        for (float j = blb_ind.y; j <= trt_ind.y; j++){
            for (float k = blb_ind.z; k <= trt_ind.z; k++){
                sp_map[glm::vec3(i, j, k)].erase(e);
            }
        }
    }

}

// void SpatialMap::update(std::shared_ptr<Boid> b){

// }



