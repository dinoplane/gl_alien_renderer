#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/random.hpp>

#include <shader_s.hpp>
#include <cube.hpp>
#include <pyramid.hpp>
#include <plane.hpp>
#include <boid.hpp>
#include <camera.hpp>

#include <spatial_map.hpp>


#include <cstdlib>
#include <thread>
#include <barrier>
#include <chrono>
#include <cassert>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

const unsigned int NUM_BOIDS = 100;

const unsigned int NUM_THREADS = 16;
const int CHUNK_SIZE = std::max((NUM_BOIDS + (NUM_THREADS - 1)) / NUM_THREADS, (unsigned int) 1);
std::barrier sync_point(NUM_THREADS);

// timing
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;

bool firstMouse = true;

float yaw   = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
float pitch =  0.0f;
float lastX =  800.0f / 2.0;
float lastY =  600.0 / 2.0;

float lightSpeed = 10.0;

Camera camera((float) SCR_WIDTH, (float) SCR_HEIGHT);


glm::vec3 lightPos = glm::vec3(0.0, 15.0, 5.0);

int setupGLFW(GLFWwindow* &window){
	// glfw: initialize and configure
    // ------------------------------
    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	#endif

    // glfw window creation
    // --------------------
    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, "LearnOpenGL", NULL, NULL);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    glfwMakeContextCurrent(window);
    glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	return 0;
}

int setupGLAD(){
	// glad: load all OpenGL function pointers
    // ---------------------------------------
    if (!gladLoadGLLoader((GLADloadproc)glfwGetProcAddress))
    {
        std::cout << "Failed to initialize GLAD" << std::endl;
        return -1;
    }

	return 0;
}

void setupShaders(unsigned int &shaderProgram){

}

// Assume there are alpha channels... need to find way to support jpeg later
std::vector<unsigned int> loadAllTextures(std::convertible_to<std::string_view> auto&& ...s){
    std::initializer_list<std::string_view> filenames = {s... };

    std::vector<unsigned int> textureNums;
    int width, height, nrChannels;
    unsigned char *data;

    size_t i = 0;
    for (auto f_n : filenames){
        std::cout <<  f_n << " ";
        unsigned int texture;
        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);
        // set the texture wrapping parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
        // set texture filtering parameters
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

        stbi_set_flip_vertically_on_load(true); // tell stb_image.h to flip loaded texture's on the y-axis.

        // load image, create texture and generate mipmaps
        data = stbi_load(f_n.data(), &width, &height, &nrChannels, 0);
        if (data)
        {
            glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
            glGenerateMipmap(GL_TEXTURE_2D);
        }
        else
        {
            std::cout << "Failed to load texture" << std::endl;
        }
        stbi_image_free(data);
        textureNums.push_back(texture);
        i += 1;
    }

    std::cout << std::endl;

    return textureNums;
}

void mouse_callback(GLFWwindow* window, double xpos, double ypos)
{
    if (glfwGetInputMode(window, GLFW_CURSOR) != GLFW_CURSOR_DISABLED){
        return;
    }
    if (firstMouse)
    {
        lastX = xpos;
        lastY = ypos;
        firstMouse = false;
    }

    float xoffset = xpos - lastX;
    float yoffset = lastY - ypos;
    lastX = xpos;
    lastY = ypos;

    camera.processMouseMovement(xoffset, yoffset);
}


void updateBoids(std::vector<std::shared_ptr<Boid>> *boids,
        unsigned int start, unsigned int end, unsigned int tid){
    for (int i = start; i < end; i++){
        (*boids)[i]->calculateForce();
        // std::string msg = std::to_string(tid) + "  Completed force!\n";
    }
        // std::cout << msg;
    sync_point.arrive_and_wait();
    for (int i = start; i < end; i++){
            (*boids)[i]->updatePosition();
            (*boids)[i]->updateMapEntry();
        // msg = std::to_string(tid) + "  Completed position!\n";
        // std::cout << msg;

    }
    sync_point.arrive_and_wait();
}

void updateBoidsForce(std::vector<std::shared_ptr<Boid>> *boids,
        unsigned int start, unsigned int end, unsigned int tid){
    for (int i = start; i < end; i++){
        if (i < NUM_BOIDS)
            (*boids)[i]->calculateForce();
        // std::string msg = std::to_string(tid) + "  Completed force!\n";
        // std::cout << msg;
    }
    sync_point.arrive_and_wait();
}

void updateBoidsPos(std::vector<std::shared_ptr<Boid>> *boids,
        unsigned int start, unsigned int end, unsigned int tid){
    for (int i = start; i < end; i++){
        if (i < NUM_BOIDS)
           (*boids)[i]->updatePosition();
        // msg = std::to_string(tid) + "  Completed position!\n";
        // std::cout << msg;
    }
    sync_point.arrive_and_wait();
}

void updateBoidsEntry(std::vector<std::shared_ptr<Boid>> *boids,
        unsigned int start, unsigned int end, unsigned int tid){
    for (int i = start; i < end; i++){
        if (i < NUM_BOIDS)
           (*boids)[i]->updateMapEntry();
        // msg = std::to_string(tid) + "  Completed position!\n";
        // std::cout << msg;
    }
    sync_point.arrive_and_wait();
}


// Add load texture function
int main()
{


	GLFWwindow* window;

    // i = std::make_shared<Boid>();


	if (setupGLFW(window) < 0){
		return -1;
	}

	if (setupGLAD() < 0){
		return -1;
	}

    // glEnable(GL_DEPTH_TEST);

    // testMap();
    // assert(false);


    // build and compile our shader program
    // ------------------------------------
    Shader ourShader("./shader/exv.vs", "./shader/exf.fs");
    Shader boidShader("./shader/boid.vs", "./shader/boid.fs");
    Shader seaShader("./shader/sea.vs", "./shader/sea.fs");



    // uncomment this call to draw in wireframe polygons.
//    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    std::vector<unsigned int> textureID = loadAllTextures("sky.png", "awesomeface.png", "noise.png");

    // tell opengl for each sampler to which texture unit it belongs to (only has to be done once)
    // -------------------------------------------------------------------------------------------
    ourShader.use();
    ourShader.setInt("texture1", 0);
    ourShader.setInt("texture2", 1);

    boidShader.use();
    boidShader.setInt("texture1", 0);
    boidShader.setInt("texture2", 1);
    boidShader.setInt("noise1", 2);

    seaShader.use();
    seaShader.setInt("texture1", 0);
    seaShader.setInt("texture2", 1);
    // pass projection matrix to shader (as projection matrix rarely changes there's no need to do this per frame)
    // -----------------------------------------------------------------------------------------------------------
    // ourShader.use();
    // glm::mat4 projection    = glm::mat4(1.0f);


    // projection = glm::perspective(glm::radians(60.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 2000.0f);
    glm::mat4 view = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 projection = camera.getProjMatrix();

    // // camera/view transformation
    Cube test;
    test.init();
    test.scale(glm::vec3(100.0));

    Cube xcoord;
    xcoord.init();
    xcoord.translate(glm::vec3(10.0, 0.0, 0.0));

    Cube ycoord;
    ycoord.init();
    ycoord.translate(glm::vec3(0.0, 10.0, 0.0));

    Cube zcoord;
    zcoord.init();
    zcoord.translate(glm::vec3(0.0, 0.0, 10.0));


    Cube light;
    light.init();
    lightPos = glm::vec3(0.0, 15.0, 5.0);
    light.translate(lightPos);

    Plane sea;
    sea.setTiles(glm::ivec2(128));
    sea.init();
    sea.scale(glm::vec3(50.0));

    // sea.printVertices();
    // assert(false);


    // AHAHAHA Smart pointers
    std::vector<std::shared_ptr<Boid>> boids;

    for (int i = 0; i < NUM_BOIDS; i++){
        boids.push_back( std::make_shared<Boid>(
                            glm::linearRand(
                                glm::vec3(-7.1f, 0.0f, -7.1f),
                                glm::vec3(7.1, 0.0f, 7.1)),
                            0.15f*glm::normalize(glm::linearRand(
                                glm::vec3(-0.15f, -0.15f, -0.15f),
                                glm::vec3(0.15f, 0.15f, 0.15f)))));
    }
    // for (auto p = Boid::boids.begin(); p != Boid::boids.end(); p++){
    //     std::cout  << "ID" << (*p)->ID  << ": "<< glm::to_string((*p)->position) << std::endl;
    // }

    // std::cout << glm::to_string(b.model.modelMat) << std::endl;
    glEnable(GL_DEPTH_TEST);
    glfwSetCursorPosCallback(window, mouse_callback);


    // timing variables
    auto begin = std::chrono::high_resolution_clock::now();
    auto start = std::chrono::high_resolution_clock::now();
    auto chkpt = std::chrono::high_resolution_clock::now();
    auto stop = std::chrono::high_resolution_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);
    std::cout << "Time taken by function: "
         << duration.count() << " microseconds" << std::endl;
    double counter = 0;
    unsigned int accum = 0;
    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {
        start = std::chrono::high_resolution_clock::now();
        // input
        // -----
        processInput(window);
        light.reset();
        light.translate(lightPos);
        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // bind textures on corresponding texture units
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureID[0]);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureID[1]);
        glActiveTexture(GL_TEXTURE2);
        glBindTexture(GL_TEXTURE_2D, textureID[2]);

        // Set up projection and view matrix
        projection = camera.getProjMatrix();
        view = camera.getViewMatrix();


        // Setup markers and area

        ourShader.use();
        ourShader.setMat4("projection", projection);
        ourShader.setMat4("view", view);

        ourShader.setInt("selected", 1);
        ourShader.setMat4("model", test.modelMat);
        // test.render();

        // ourShader.use();
        // ourShader.setMat4("model", xcoord.modelMat);
        // xcoord.render();

        // ourShader.use();
        // ourShader.setMat4("model", ycoord.modelMat);
        // ycoord.render();

        // ourShader.use();
        // ourShader.setMat4("model", zcoord.modelMat);
        // zcoord.render();

        ourShader.use();
        ourShader.setMat4("model", light.modelMat);
        light.render();

        seaShader.use();
        seaShader.setMat4("projection", projection);
        seaShader.setMat4("view", view);
        seaShader.setInt("selected", 0);
        seaShader.setFloat("uTime", glfwGetTime());
        seaShader.setMat4("model", sea.modelMat);
        seaShader.setVec3("uLightPos", lightPos);
        seaShader.setVec3("uViewPos", camera.position);


        sea.render();


        std::thread *thread_array = new std::thread[NUM_THREADS];
        /*
        { // Force Calculation
            chkpt = std::chrono::high_resolution_clock::now();
            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                unsigned int start = CHUNK_SIZE * tid;
                // unsigned int end = start + CHUNK_SIZE;
                unsigned int end = std::min(start + CHUNK_SIZE, NUM_BOIDS);

                // if (start < end)
                // std::cout << start << " "  << end << std::endl;
                thread_array[tid] =
                    std::thread(updateBoidsForce, &boids, start, end, tid);
            }

            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                thread_array[tid].join();
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<
                        std::chrono::microseconds>(stop - chkpt);

            std::cout << "Time taken by Update Force: "
                << duration.count() << " microseconds" << std::endl;
        }

        { // Position Update
            chkpt = std::chrono::high_resolution_clock::now();
            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                unsigned int start = CHUNK_SIZE * tid;
                // unsigned int end = start + CHUNK_SIZE;
                unsigned int end = std::min(start + CHUNK_SIZE, NUM_BOIDS);

                // if (start < end)
                // std::cout << start << " "  << end << std::endl;
                thread_array[tid] =
                    std::thread(updateBoidsPos, &boids, start, end, tid);
            }

            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                thread_array[tid].join();
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<
                        std::chrono::microseconds>(stop - chkpt);

            std::cout << "Time taken by Update Position: "
                << duration.count() << " microseconds" << std::endl;
        }

        { // Map Entry Update
            chkpt = std::chrono::high_resolution_clock::now();
            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                unsigned int start = CHUNK_SIZE * tid;
                // unsigned int end = start + CHUNK_SIZE;
                unsigned int end = std::min(start + CHUNK_SIZE, NUM_BOIDS);

                // if (start < end)
                // std::cout << start << " "  << end << std::endl;
                thread_array[tid] =
                    std::thread(updateBoidsEntry, &boids, start, end, tid);
            }

            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                thread_array[tid].join();
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<
                        std::chrono::microseconds>(stop - chkpt);

            std::cout << "Time taken by Update Entry: "
                << duration.count() << " microseconds" << std::endl;
        }

        // for (int i = 0; i < NUM_BOIDS; i++){
        //     boids[i]->calculateForce();
        // }
        */

        { // All calculation
            chkpt = std::chrono::high_resolution_clock::now();
            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                unsigned int start = CHUNK_SIZE * tid;
                // unsigned int end = start + CHUNK_SIZE;
                unsigned int end = std::min(start + CHUNK_SIZE, NUM_BOIDS);

                // if (start < end)
                // std::cout << start << " "  << end << std::endl;
                thread_array[tid] =
                    std::thread(updateBoids, &boids, start, end, tid);
            }

            for (unsigned int tid = 0; tid < NUM_THREADS; tid++){
                thread_array[tid].join();
            }

            stop = std::chrono::high_resolution_clock::now();
            duration = std::chrono::duration_cast<
                        std::chrono::microseconds>(stop - chkpt);

            std::cout << "Time taken by Update Boids: "
                << duration.count() << " microseconds" << std::endl;
        }

        { // Rendering
            chkpt = std::chrono::high_resolution_clock::now();
            boidShader.use();
            boidShader.setMat4("projection", projection);
            boidShader.setInt("selected", 1);
            boidShader.setMat4("view", view);

            for (int i = 0; i < NUM_BOIDS; i++){
                // boids[i]->updatePosition();
                boidShader.use();
                boidShader.setMat4("model", boids[i]->model.modelMat);
                boids[i]->render();
            }

            // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
            // -------------------------------------------------------------------------------
            glfwSwapBuffers(window);
            glfwPollEvents();
            deltaTime = glfwGetTime() - lastFrame;
            lastFrame = glfwGetTime();

            stop = std::chrono::high_resolution_clock::now();

            duration = std::chrono::duration_cast<
                        std::chrono::microseconds>(stop - chkpt);

            std::cout << "Time taken by Rendering: "
                << duration.count() << " microseconds" << std::endl;
        }


        duration = std::chrono::duration_cast<std::chrono::microseconds>(stop - start);

        std::cout << "Time taken by function: "
            << duration.count() << " microseconds" << std::endl;
        counter++;
        accum += duration.count();
    }

    std::cout << "Average Time taken by function: "
        << (accum/counter) << " microseconds" << std::endl;


    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    test.deleteBuffers();
    ourShader.deleteProgram();
    boidShader.deleteProgram();
    seaShader.deleteProgram();


    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();
    return 0;
}


// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------

float SPEED_MULTIPLIER = 10;
void processInput(GLFWwindow *window)
{
    double mWidth = (double) SCR_WIDTH;
    double mHeight = (double) SCR_HEIGHT;

    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_NORMAL){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }

    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }


    // std::cout << mWidth << " " << mHeight << std::endl;
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS){
        Boid::boid_map.printMap();
        glfwSetWindowShouldClose(window, true);

    }

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
        // std::cout << "PRESSED " << deltaTime << std::endl;
        camera.moveForward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
        camera.moveBackward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
        camera.moveLeft(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
        camera.moveRight(deltaTime);
    }

    if (glfwGetKey(window, GLFW_KEY_UP) == GLFW_PRESS){
        // std::cout << "PRESSED " << deltaTime << std::endl;
        lightPos += glm::vec3(0.0, 0.0, -1.0) * lightSpeed * deltaTime;
    }
    if (glfwGetKey(window, GLFW_KEY_DOWN) == GLFW_PRESS){
        lightPos += glm::vec3(0.0, 0.0, 1.0) * lightSpeed * deltaTime;
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT) == GLFW_PRESS){
        lightPos += glm::vec3(-1.0, 0.0, 0.0) * lightSpeed * deltaTime;
    }
    if (glfwGetKey(window, GLFW_KEY_RIGHT) == GLFW_PRESS){
        lightPos += glm::vec3(1.0, 0.0, 0.0) * lightSpeed * deltaTime;
    }
    if (glfwGetKey(window, GLFW_KEY_SLASH) == GLFW_PRESS){
        lightPos += glm::vec3(0.0, 1.0, 0.0) * lightSpeed * deltaTime;
    }
    if (glfwGetKey(window, GLFW_KEY_PERIOD) == GLFW_PRESS){
        lightPos += glm::vec3(0.0, -1.0, 0.0) * lightSpeed * deltaTime;
    }



}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
    camera.setPerspectiveSize(width, height);
}

