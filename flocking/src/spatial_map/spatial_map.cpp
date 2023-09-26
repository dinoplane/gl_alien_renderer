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

// Inspiration from: https://www.youtube.com/watch?v=sx4IIQL0x7c

/*
    Returns the index of the cell the position is in.

*/
glm::ivec3 SpatialMap::_key(glm::vec3 pos){
    return  glm::sign(pos)*glm::floor(glm::abs(pos /(unit_dims)));
}

void SpatialMap::printMap(){
    printMap(minbound, maxbound);

    // std::cout << "Out of Bounds: " << sp_map[glm::vec3(-1.0f)].size() << std::endl;
}

void SpatialMap::printMap(glm::ivec3 blb_ind, glm::ivec3 trt_ind){
    std::cout << "grid for " << glm::to_string(blb_ind) << " to " << glm::to_string(trt_ind) << std::endl;
    size_t overall_count = 0;
    for (int j = blb_ind.y; j <= trt_ind.y; j++){
        std::cout << "Y = " << j << std::endl;
        std::cout << "     ";
        size_t layer_count = 0;
        for (int i = blb_ind.x; i <= trt_ind.x; i++){
            std::printf(" %03d ", i);
        }
        std::cout << std::endl;
        for (int k = blb_ind.z; k <= trt_ind.z; k++){
            std::printf("%03d  ", k);
            for (int i = blb_ind.x; i <= trt_ind.x; i++){
                glm::ivec3 key = glm::ivec3(i, j, k);
                std::printf(" %03lu ", sp_map[key].s_set.size());
                // std::unordered_map<glm::vec3, int> sp_set;
                // sp_set[key] = 1;
                layer_count += sp_map[key].s_set.size();
                // std::cout << glm::to_string(key) << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "Total in layer: " << layer_count << std::endl;
        std::cout << std::endl;
        overall_count += layer_count;
    }
    std::cout << "Overall Count: " << overall_count << std::endl;

    // std::cout << "Out of Bounds: " << sp_map[glm::vec3(-1.0f)].size() << std::endl;
}

void SpatialMap::displayContents(){
    displayContents(minbound, maxbound);

    // std::cout << "Out of Bounds: " << sp_map[glm::vec3(-1.0f)].size() << std::endl;
}

void SpatialMap::displayContents(glm::ivec3 blb_ind, glm::ivec3 trt_ind){

    std::cout << "------------------------------------------------------------------------" << std::endl;
    std::cout << "------------------------------------------------------------------------" << std::endl;

    std::cout << "contents for " << glm::to_string(blb_ind) << " to " << glm::to_string(trt_ind) << std::endl;

    for (int k = blb_ind.z; k <= trt_ind.z; k++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int i = blb_ind.x; i <= trt_ind.x; i++){
                glm::ivec3 key = glm::ivec3(i, j, k);

                std::cout << "---------------------------------------------" << std::endl;
                std::cout << glm::to_string(key) << std::endl;
                std::cout << "---------------------------------------------" << std::endl;

                for (std::shared_ptr<SpatialEntry> e : sp_map[key].s_set){
                    std::cout << "entry pos: " << glm::to_string(e->position) << std::endl;
                    std::cout << "dimension: " << glm::to_string(e->dimensions) << std::endl;
                    std::cout << "BLB indic: " << glm::to_string(e->blb_ind) << std::endl;
                    std::cout << "TRT indic: " << glm::to_string(e->trt_ind) << std::endl;
                }

                // std::cout << glm::to_string(key) << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
}


std::shared_ptr<SpatialEntry> SpatialMap::insert(Boid* b){
    // Everything below can be private?
    auto entry = std::make_shared<SpatialEntry>(b, b->position, b->dimensions);

    return _insert(entry);
}

std::shared_ptr<SpatialEntry> SpatialMap::_insert(std::shared_ptr<SpatialEntry> e){
    glm::vec3 blb = e->position - e->dimensions/2.0f;
    glm::vec3 trt = e->position + e->dimensions/2.0f;

    glm::ivec3 blb_ind = _key(blb);
    glm::ivec3 trt_ind = _key(trt);

    e->blb_ind = blb_ind;
    e->trt_ind = trt_ind;

    // std::cout  << "Bottom Left Back: " << glm::to_string(blb_ind) << std::endl;
    // std::cout << "Top Right Top: " << glm::to_string(trt_ind) << std::endl;

    /*
    map_lock.lock();
    bool blb_out = glm::any(glm::lessThanEqual(blb_ind, minbound));
    bool trt_out =glm::any(glm::greaterThanEqual(trt_ind, maxbound));
    glm::ivec3 newbound;
    if (blb_out){
        newbound = glm::ivec3(glm::lessThanEqual(blb_ind, minbound)) * blb_ind +
            glm::ivec3(glm::greaterThanEqual(blb_ind, minbound)) * minbound;
        for (int i = newbound.x; i < minbound.x; i++){
            for (int j = newbound.y; j < minbound.y; j++){
                for (int k = newbound.z; k < minbound.z; k++){
                    glm::ivec3 key = glm::ivec3(i, j, k);
                    sp_map[key];
                }
            }
        }
        minbound = newbound;
    }

    if (trt_out){
        newbound = glm::ivec3(glm::greaterThanEqual(trt_ind, maxbound)) * trt_ind +
            glm::ivec3(glm::lessThan(trt_ind, maxbound)) * maxbound;

        for (int i = maxbound.x + 1; i <= newbound.x; i++){
            for (int j = maxbound.y + 1; j <= newbound.y; j++){
                for (int k = maxbound.z + 1; k <= newbound.z; k++){
                    glm::ivec3 key = glm::ivec3(i, j, k);
                    sp_map[key];
                }
            }
        }
        maxbound = newbound;
    }
    map_lock.unlock();
    */

    for (int i = blb_ind.x; i <= trt_ind.x; i++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int k = blb_ind.z; k <= trt_ind.z; k++){
                glm::ivec3 key = glm::ivec3(i, j, k);

                map_lock.lock(); // Naive solution...
                // sp_map[key].set_lock.lock();
                sp_map[key].s_set.emplace(e);
                // sp_map[key].set_lock.unlock();
                map_lock.unlock();


            }
        }
    }

    return e;
}


void SpatialMap::remove(std::shared_ptr<SpatialEntry> e){
    glm::ivec3 blb_ind = e->blb_ind;
    glm::ivec3 trt_ind = e->trt_ind;


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


    for (int i = blb_ind.x; i <= trt_ind.x; i++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int k = blb_ind.z; k <= trt_ind.z; k++){
                glm::ivec3 key = glm::ivec3(i, j, k);

                map_lock.lock();
                // sp_map[key].set_lock.lock();
                sp_map[key].s_set.erase(e);
                // sp_map[key].set_lock.unlock();
                map_lock.unlock();
            }
        }
    }

}

void SpatialMap::update(std::shared_ptr<SpatialEntry>e){
    glm::vec3 blb = e->boid->position - e->boid->dimensions/2.0f;
    glm::vec3 trt = e->boid->position + e->boid->dimensions/2.0f;

    glm::ivec3 blb_ind = _key(blb);
    glm::ivec3 trt_ind = _key(trt);

    if (blb_ind != e->blb_ind || trt_ind != e->trt_ind){
        remove(e);
        e->position = e->boid->position;
        _insert(e);
    } else {
        e->position = e->boid->position;
    }


}

std::unordered_set<Boid*> SpatialMap::getNearby(Boid* b, float range){
    glm::vec3 src = b->position;
    std::unordered_set<Boid*> ret;
    glm::vec3 blb = src - range;
    glm::vec3 trt = src + range;

    glm::ivec3 blb_ind = _key(blb);
    glm::ivec3 trt_ind = _key(trt);


    // std::cout  << "Bottom Left Back: " << glm::to_string(blb_ind) << std::endl;
    // std::cout << "Top Right Top: " << glm::to_string(trt_ind) << std::endl;

    for (int i = blb_ind.x; i <= trt_ind.x; i++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int k = blb_ind.z; k <= trt_ind.z; k++){
                glm::vec3 key = glm::vec3(i, j, k);

                if (sp_map.contains(key)){
                    auto bgn = sp_map[key].s_set.begin();
                    auto end = sp_map[key].s_set.end();
                    for (; bgn != end; bgn++){
                        // std::cout << glm::to_string((*bgn)->position) <<
                        // " vs " <<
                        //  glm::to_string((*bgn)->boid->position) << std::endl ;
                        float dist = glm::length((*bgn)->position - src);
                        if (b->ID != (*bgn)->boid->ID &&
                            0.0f < dist &&
                            dist <= range)
                            ret.emplace((*bgn)->boid);
                    }
                }
            }
        }
    }
    return ret;
}

// std::unordered_set<std::shared_ptr<SpatialEntry>> getNearby(glm::vec3 src, float range){
//     std::unordered_set<std::shared_ptr<SpatialEntry>> ret;

