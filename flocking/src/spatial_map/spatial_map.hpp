#ifndef SPATIAL_MAP
#define SPATIAL_MAP


#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>

#define GLM_ENABLE_EXPERIMENTAL
#include <glm/gtx/hash.hpp>

#include <unordered_map>
#include <vector>
#include <unordered_set> // consider and benchmark unordered
#include <memory>

#include <boid.hpp>


#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

struct SpatialEntry {
    std::shared_ptr<Boid> boid;
    glm::vec3 position;
    glm::vec3 dimensions;
    glm::vec3 blb_ind;
    glm::vec3 trt_int;
    // unsigned long hash(){return std::hash<glm::vec3>(position);};
};

class SpatialMap {
    public:
    // Let's generalize this later
        std::unordered_map<glm::vec3, std::unordered_set<std::shared_ptr<SpatialEntry>>> sp_map;

        glm::vec3 _key(glm::vec3 pos);
        void _insert(std::shared_ptr<SpatialEntry> e);
        // void _remove(std::shared_ptr<Boid> b);
        // void _update(std::shared_ptr<Boid> b);// Update only if needed!

    public:
        glm::vec3 bound_dims; // the total size
        glm::vec3 minbound;
        glm::vec3 maxbound;
        glm::vec3 unit_dims;
        glm::vec3 grid_dims; // number of cubes along each axis

    SpatialMap(glm::vec3 dims=glm::vec3(10.0f), glm::vec3 gsize=glm::vec3(3.0f)){ // We could test the affects of granulity
        bound_dims = dims;
        grid_dims = gsize;

        glm::vec3 blb = glm::vec3(-dims.x/2, -dims.y/2, -dims.z/2);
        glm::vec3 udims = glm::vec3(
                                    dims.x/gsize.x,
                                    dims.y/gsize.y,
                                    dims.z/gsize.z
                                );
        for (int i = 0; i < gsize.x; i++){
            for (int j = 0; j < gsize.y; j++){
                for (int k = 0; k < gsize.z; k++){

                    glm::vec3 key = glm::vec3(i, j, k);

                    sp_map[key] = std::unordered_set<std::shared_ptr<SpatialEntry>>();
                    // std::unordered_map<glm::vec3, int> sp_set;
                    // sp_set[key] = 1;

                    // std::cout << glm::to_string(key) << std::endl;
                }
            }
        }
        sp_map[glm::vec3(-1.0f)] = std::unordered_set<std::shared_ptr<SpatialEntry>>();

        minbound = blb;
        maxbound = blb + gsize*udims;
        unit_dims = udims;

        // std::cout << "Minbound: " << glm::to_string(minbound) << std::endl;
        // std::cout << "Maxbound: " << glm::to_string(maxbound) << std::endl;
    }

    ~SpatialMap(){

    }


    SpatialMap(const SpatialMap& other){
        bound_dims = other.bound_dims;
        grid_dims = other.grid_dims;

    }
    SpatialMap(SpatialMap&& other){
        bound_dims = other.bound_dims;
        grid_dims = other.grid_dims;
    }

    SpatialMap& operator=(const SpatialMap& other){
        bound_dims = other.bound_dims;
        grid_dims = other.grid_dims;
        return *this;
    }
    SpatialMap& operator=(SpatialMap&& other){
        bound_dims = std::move(other.bound_dims);
        grid_dims = std::move(other.grid_dims);
        return *this;
    }


    void printMap();
    void insert(Boid* b);
    void remove(std::shared_ptr<SpatialEntry> e);
    // void update(std::shared_ptr<Boid> b);// Update only if needed!
    // std::unordered_set<std::shared_ptr<Boid>> getNear(7dfloat range);
};

#endif