#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/random.hpp>

#include <shader_s.hpp>
#include <camera.hpp>
#include <mesh.hpp>

#include <scene.hpp>
#include <renderer.hpp>

#include <cstdlib>
#include <thread>
#include <barrier>
#include <chrono>
#include <cassert>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>
// #include <fmt/core.h>
// #include <fmt/printf.h>



void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window, Camera *camera);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;

// int CHUNK_SIZE = std::max((NUM_BOIDS + (NUM_THREADS - 1)) / NUM_THREADS, (unsigned int) 1);
// // std::barrier sync_point(NUM_THREADS);

// // timing
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;

bool firstMouse = true;

// float yaw   = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
// float pitch =  0.0f;
float lastX =  800.0f / 2.0;
float lastY =  600.0 / 2.0;

// float lightSpeed = 10.0;
Scene scene;

// glm::vec3 lightPos = glm::vec3(0.0, 15.0, 5.0);

void setupCmdArgs(int argc, char **argv){
    // int i = 1;
    // while (i < argc){
    //     if (argv[i][0] != '-') {
    //         std::cout << "Usage: flocking -b NUM_BOUNDS -t NUM_THREADS\n";
    //         exit(EXIT_FAILURE);
    //     }
    //     switch (argv[i][1]){
    //         case 'b':
    //             i++;
    //             NUM_BOIDS = std::stoi(argv[i]);
    //             break;
    //         case 't':
    //             i++;
    //             NUM_THREADS = std::stoi(argv[i]);
    //             break;
    //     }
    //     i++;
    // }
    // CHUNK_SIZE = std::max((NUM_BOIDS + (NUM_THREADS - 1)) / NUM_THREADS, (unsigned int) 1);

}

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

unsigned int loadCubeMap(std::convertible_to<std::string_view> auto&& ...s){
    std::initializer_list<std::string_view> textures_faces = {s... };

    unsigned int textureID;
    int width, height, nrChannels;
    unsigned char *data;

    glGenTextures(1, &textureID);
    glBindTexture(GL_TEXTURE_CUBE_MAP, textureID);
    int i = 0;
    for(auto f_n : textures_faces)
    {
        data = stbi_load(f_n.data(), &width, &height, &nrChannels, 0);

        if (data){
            glTexImage2D(
                GL_TEXTURE_CUBE_MAP_POSITIVE_X + i,
                0, GL_RGB, width, height, 0, GL_RGB, GL_UNSIGNED_BYTE, data
            );
        } else {
            std::cout << "Failed to load texture "  << f_n.data()  << std::endl;
        }
        stbi_image_free(data);
        i++;
    }

    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_CUBE_MAP, GL_TEXTURE_WRAP_R, GL_CLAMP_TO_EDGE);

    return textureID;
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

        scene.cameras[0].processMouseMovement(xoffset, yoffset);
}
template<typename T>
void printVector(const std::vector<T>& vec) {
    std::cout << "Elements of the vector:" << std::endl;
    for (const T& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

// Add load texture function
int main(int argc, char **argv)
{

    setupCmdArgs(argc,argv);


	GLFWwindow* window;

    // i = std::make_shared<Boid>();


	if (setupGLFW(window) < 0){
		return -1;
	}

	if (setupGLAD() < 0){
		return -1;
	}


    // build and compile our shader program
    // ------------------------------------


    // uncomment this call to draw in wireframe polygons.
//    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    // std::vector<unsigned int> textureID = loadAllTextures("./resources/assets/sky.png", "./resources/assets/awesomeface.png", "./resources/assets/noise.png");
    // unsigned int skyBoxID = loadCubeMap(
    //                         "./resources/assets/bluecloud_rt.jpg",
    //                         "./resources/assets/bluecloud_lf.jpg",
    //                         "./resources/assets/bluecloud_up.jpg",
    //                         "./resources/assets/bluecloud_dn.jpg",
    //                         "./resources/assets/bluecloud_bk.jpg",
    //                         "./resources/assets/bluecloud_ft.jpg"
    //                     );
    // tell opengl for each sampler to which texture unit it belongs to (only has to be done once)
    // -------------------------------------------------------------------------------------------

    scene = Scene::GenerateDefaultScene();

    // printVector(scene.meshes);
    // fmt::printf("Scene generated %u\n", scene.meshes.size());


    // assert(false);
    Renderer renderer = Renderer();

    // projection = glm::perspective(glm::radians(60.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 2000.0f);

    // sea.printVertices();
    // assert(false);
    glEnable(GL_DEPTH_TEST);

    // auto mouseCallback = [&](GLFWwindow* window, double xpos, double ypos) -> void {

    // };

    glfwSetCursorPosCallback(window, mouse_callback);


    // std::cout << glm::to_string(b.model.modelMat) << std::endl;
    // glEnable(GL_DEPTH_TEST);
    // glfwSetCursorPosCallback(window, mouse_callback);

    /*

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
    */
    while (!glfwWindowShouldClose(window))
    {
        // input
        // -----
        processInput(window, &scene.cameras[0]);

        renderer.Render(&scene);

        // light.reset();
        // light.translate(lightPos);
        // render
        // ------
        // glClearColor(0.0f, 0.0f, 0.0f, 1.0f);
        // glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // // bind textures on corresponding texture units
        // glActiveTexture(GL_TEXTURE0);
        // glBindTexture(GL_TEXTURE_2D, textureID[0]);
        // glActiveTexture(GL_TEXTURE1);
        // glBindTexture(GL_TEXTURE_2D, textureID[1]);
        // glActiveTexture(GL_TEXTURE2);
        // glBindTexture(GL_TEXTURE_2D, textureID[2]);
        // glDepthMask(GL_FALSE);
        // glBindTexture(GL_TEXTURE_CUBE_MAP, skyBoxID);


        // Set up projection and view matrix
        // projection = camera.getProjMatrix();
        // view = camera.getViewMatrix();

        // skyboxShader.use();
        // skyboxShader.setMat4("projection", projection);
        // skyboxShader.setMat4("view", view);
        // skybox.render();

        glDepthMask(GL_TRUE);

        // Setup markers and area



        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
        deltaTime = glfwGetTime() - lastFrame;
        lastFrame = glfwGetTime();

    }



    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    // test.deleteBuffers();
    // ourShader.deleteProgram();
    // boidShader.deleteProgram();
    // // baseShader.deleteProgram();
    // seaShader.deleteProgram();



    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();

   return 0;
}


// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------


void framebuffer_size_callback(GLFWwindow* window, int width, int height){
    //
}

float SPEED_MULTIPLIER = 10;
void processInput(GLFWwindow *window, Camera *camera){
    double mWidth = (double) SCR_WIDTH;
    double mHeight = (double) SCR_HEIGHT;

    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_NORMAL){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }

    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }

    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS){
        glfwSetWindowShouldClose(window, true);

    }

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
        std::cout << "PRESSED " << deltaTime << std::endl;
        camera->moveForward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
        camera->moveBackward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
        camera->moveLeft(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
        camera->moveRight(deltaTime);
    }
}


/*
{

    // std::cout << mWidth << " " << mHeight << std::endl;


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
*/

