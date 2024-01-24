#ifndef SPATIAL_MAP
#define SPATIAL_MAP


#include <glm/glm.hpp>
#include <glm/gtx/norm.hpp>
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

#include <mutex>

class Boid;

struct SpatialEntry {
    Boid* boid;
    glm::vec3 position;
    glm::vec3 dimensions;
    glm::ivec3 blb_ind;
    glm::ivec3 trt_ind;
    // unsigned long hash(){return std::hash<glm::vec3>(position);};
};


// Basically a wrapper of a set with a mutex
struct SpatialSet{
    std::mutex set_lock;
    std::unordered_set<std::shared_ptr<SpatialEntry>> s_set;
};

class SpatialMap {
    public:
        std::mutex map_lock;
    // Let's generalize this later
        std::unordered_map<glm::ivec3, SpatialSet> spatial_map;

        glm::ivec3 _key(glm::vec3 pos);
        std::shared_ptr<SpatialEntry> _insert(std::shared_ptr<SpatialEntry> e);
        // void _remove(Boid* b);
        // void _update(Boid* b);// Update only if needed!

    public:
        // glm::vec3 bound_dims; // the total size
        glm::ivec3 minbound;
        glm::ivec3 maxbound;
        glm::vec3 unit_dims;
        glm::ivec3 grid_dims; // number of cubes along each axis

    // SpatialMap(glm::vec3 dims=glm::vec3(20.0f), glm::vec3 gsize=glm::vec3(3.0f)){ // We could test the affects of granulity
    //     // bound_dims = dims;
    //     grid_dims = gsize;

    //     glm::vec3 blb = glm::vec3(-dims.x/2, -dims.y/2, -dims.z/2);
    //     glm::vec3 udims = glm::vec3(
    //                                 dims.x/gsize.x,
    //                                 dims.y/gsize.y,
    //                                 dims.z/gsize.z
    //                             );

    //     maxbound = glm::vec3((int) gsize.x/2, (int) gsize.y/2, (int) gsize.z/2);
    //     minbound = -maxbound;
    //     for (int i = -gsize.x/2; i < gsize.x/2; i++){
    //         for (int j = -gsize.y/2; j < gsize.y/2; j++){
    //             for (int k = -gsize.z/2; k < gsize.z/2; k++){

    //                 glm::vec3 key = glm::vec3(i, j, k);

    //                 spatial_map[key].s_set = std::unordered_set<std::shared_ptr<SpatialEntry>>();
    //                 // std::unordered_map<glm::vec3, int> spatial_set;
    //                 // spatial_set[key] = 1;

    //                 std::cout << glm::to_string(key) << std::endl;
    //             }
    //         }
    //     }

    //     unit_dims = udims;
    // }

    // SpatialMap(glm::vec3 dims=glm::vec3(20.0f), glm::vec3 gsize=glm::vec3(3.0f)){ // We could test the affects of granulity
    //     grid_dims = gsize;

    //     glm::vec3 blb = glm::vec3(-dims.x/2, -dims.y/2, -dims.z/2);
    //     glm::vec3 udims = glm::vec3(
    //                                 dims.x/gsize.x,
    //                                 dims.y/gsize.y,
    //                                 dims.z/gsize.z
    //                             );

    //     maxbound = glm::vec3((int) gsize.x/2, (int) gsize.y/2, (int) gsize.z/2);
    //     minbound = -maxbound;
    //     for (int i = -gsize.x/2; i < gsize.x/2; i++){
    //         for (int j = -gsize.y/2; j < gsize.y/2; j++){
    //             for (int k = -gsize.z/2; k < gsize.z/2; k++){

    //                 glm::vec3 key = glm::vec3(i, j, k);

    //                 spatial_map[key].s_set = std::unordered_set<std::shared_ptr<SpatialEntry>>();
    //                 // std::unordered_map<glm::vec3, int> spatial_set;
    //                 // spatial_set[key] = 1;

    //                 std::cout << glm::to_string(key) << std::endl;
    //             }
    //         }
    //     }

    //     unit_dims = udims;
    // }

    SpatialMap(glm::ivec3 udims=glm::ivec3(5.0f), glm::ivec3 gdims=glm::ivec3(11.0f)){ // We could test the affects of granulity
        maxbound = glm::ivec3( gdims/2);
        minbound = glm::ivec3(-gdims/2);
        unit_dims = udims;

        for (int i = -gdims.x/2; i <= gdims.x/2; i++){
            for (int j = -gdims.y/2; j <= gdims.y/2; j++){
                for (int k = -gdims.z/2; k <= gdims.z/2; k++){
                    glm::ivec3 key = glm::ivec3(i, j, k);

                    spatial_map[key];
                }
            }
        }
    }



    ~SpatialMap(){

        auto it = spatial_map.begin();

        for (;it != spatial_map.end(); ++it){
            it->second.s_set.clear();
        }
        spatial_map.clear();
    }


    SpatialMap(const SpatialMap& other){
        grid_dims = other.grid_dims;

    }
    SpatialMap(SpatialMap&& other){
        grid_dims = other.grid_dims;
    }

    SpatialMap& operator=(const SpatialMap& other){
        grid_dims = other.grid_dims;
        return *this;
    }
    SpatialMap& operator=(SpatialMap&& other){
        grid_dims = std::move(other.grid_dims);
        return *this;
    }


    void printMap();
    void printMap(glm::ivec3 blb_ind, glm::ivec3 trt_ind);

    void displayContents();
    void displayContents(glm::ivec3 blb_ind, glm::ivec3 trt_ind);


    std::shared_ptr<SpatialEntry> insert(Boid* b);
    void remove(Boid* b);
    void update(Boid* b);// Update only if needed!
    std::vector<Boid*> getNearby(Boid* b, float range);
    // SpatialSet getNearby(glm::vec3 src, float range);

};

#endif