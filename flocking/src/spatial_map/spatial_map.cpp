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
#include <atomic>


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

    // std::cout << "Out of Bounds: " << spatial_map[glm::vec3(-1.0f)].size() << std::endl;
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
                std::printf(" %03lu ", spatial_map[key].s_set.size());
                // std::unordered_map<glm::vec3, int> spatial_set;
                // spatial_set[key] = 1;
                layer_count += spatial_map[key].s_set.size();
                // std::cout << glm::to_string(key) << std::endl;
            }
            std::cout << std::endl;
        }
        std::cout << "Total in layer: " << layer_count << std::endl;
        std::cout << std::endl;
        overall_count += layer_count;
    }
    std::cout << "Overall Count: " << overall_count << std::endl;

    // std::cout << "Out of Bounds: " << spatial_map[glm::vec3(-1.0f)].size() << std::endl;
}

void SpatialMap::displayContents(){
    displayContents(minbound, maxbound);

    // std::cout << "Out of Bounds: " << spatial_map[glm::vec3(-1.0f)].size() << std::endl;
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

                for (std::shared_ptr<SpatialEntry> e : spatial_map[key].s_set){
                    std::cout << "entry pos: " << glm::to_string(e->position) << std::endl;
                    std::cout << "dimension: " << glm::to_string(e->dimensions) << std::endl;
                    std::cout << "BLB indic: " << glm::to_string(e->blb_ind) << std::endl;
                    std::cout << "TRT indic: " << glm::to_string(e->trt_ind) << std::endl;
                }
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

    for (int i = blb_ind.x; i <= trt_ind.x; i++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int k = blb_ind.z; k <= trt_ind.z; k++){
                glm::ivec3 key = glm::ivec3(i, j, k);
                while (true){
                    try {
                        const std::lock_guard<std::mutex> set_lock_guard(spatial_map.at(key).set_lock);
                        spatial_map.at(key).s_set.emplace(e);
                        break;
                    } catch (const std::out_of_range &error) {
                        const std::lock_guard<std::mutex> map_lock_guard(map_lock);
                        spatial_map[key];
                    }
                }
            }
        }
    }

    return e;
}


void SpatialMap::remove(std::shared_ptr<SpatialEntry> e){
    glm::ivec3 blb_ind = e->blb_ind;
    glm::ivec3 trt_ind = e->trt_ind;

    for (int i = blb_ind.x; i <= trt_ind.x; i++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int k = blb_ind.z; k <= trt_ind.z; k++){
                glm::ivec3 key = glm::ivec3(i, j, k);

                while (true){
                    try {
                        const std::lock_guard<std::mutex> set_lock_guard(spatial_map.at(key).set_lock);
                        spatial_map.at(key).s_set.erase(e);
                        break;
                    } catch (const std::out_of_range &error) {
                        const std::lock_guard<std::mutex> map_lock_guard(map_lock);
                        spatial_map[key];
                    }
                }
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

std::vector<Boid*> SpatialMap::getNearby(Boid* b, float range){
    glm::vec3 src = b->position;
    std::vector<Boid*> ret;
    glm::vec3 blb = src - range;
    glm::vec3 trt = src + range;

    glm::ivec3 blb_ind = _key(blb);
    glm::ivec3 trt_ind = _key(trt);

    for (int i = blb_ind.x; i <= trt_ind.x; i++){
        for (int j = blb_ind.y; j <= trt_ind.y; j++){
            for (int k = blb_ind.z; k <= trt_ind.z; k++){
                glm::vec3 key = glm::vec3(i, j, k);

                if (spatial_map.contains(key)){
                    auto bgn = spatial_map[key].s_set.begin();
                    auto end = spatial_map[key].s_set.end();
                    for (; bgn != end; bgn++){
                        float dist = glm::length2((*bgn)->position - src);
                        if (b->ID != (*bgn)->boid->ID &&
                            0.0f < dist &&
                            dist <= range)
                            ret.push_back((*bgn)->boid);
                    }
                }
            }
        }
    }
    return ret;
}

