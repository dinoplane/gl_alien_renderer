#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/random.hpp>

#include <shader_s.hpp>
#include <camera.hpp>
#include <mesh.hpp>

#include <scene.hpp>
#include <scenedata.hpp>
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




#include <stdio.h>
#include <stdlib.h>

#include "linmath.h"

static const char* vertex_shader_text =
"#version 110\n"
"uniform mat4 MVP;\n"
"attribute vec2 vPos;\n"
"varying vec2 texcoord;\n"
"void main()\n"
"{\n"
"    gl_Position = MVP * vec4(vPos, 0.0, 1.0);\n"
"    texcoord = vPos;\n"
"}\n";

static const char* fragment_shader_text =
"#version 110\n"
"uniform sampler2D texture;\n"
"uniform vec3 color;\n"
"varying vec2 texcoord;\n"
"void main()\n"
"{\n"
"    gl_FragColor = vec4(color * texture2D(texture, texcoord).rgb, 1.0);\n"
"}\n";

static const vec2 vertices[4] =
{
    { 0.f, 0.f },
    { 1.f, 0.f },
    { 1.f, 1.f },
    { 0.f, 1.f }
};




void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window, Renderer* renderer, Scene* scene);

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

SceneData sceneData;
Scene scene;
// Renderer renderer;
std::vector<Renderer> renderers;

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
int createGLFWwindow(GLFWwindow* &window, const char* windowName, GLFWwindow* share = NULL){
    // glfw window creation
    // --------------------

    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, windowName, NULL, share);
    if (window == NULL)
    {
        std::cout << "Failed to create GLFW window" << std::endl;
        glfwTerminate();
        return -1;
    }
    // glfwMakeContextCurrent(window);
    // glfwSetFramebufferSizeCallback(window, framebuffer_size_callback);
	return 0;
}

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

int setupGLFW(){
	// glfw: initialize and configure
    // ------------------------------

    glfwSetErrorCallback(error_callback);

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);

	#ifdef __APPLE__
		glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
	#endif

    return 1;
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

        scene.cameras[renderers[0].mainCameraIdx].processMouseMovement(xoffset, yoffset);
}
template<typename T>
void printVector(const std::vector<T>& vec) {
    std::cout << "Elements of the vector:" << std::endl;
    for (const T& elem : vec) {
        std::cout << elem << " ";
    }
    std::cout << std::endl;
}

void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
    if (key == GLFW_KEY_TAB && action == GLFW_PRESS){

        renderers[0].SwitchCamera(scene);
    }

    // if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
    //     glfwSetWindowShouldClose(window, true);
    // }

}


int main(int argc, char** argv)
{
    std::vector<GLFWwindow*> windows;
    windows.resize(2);
    GLuint texture, program, vertex_buffer;
    GLint mvp_location, vpos_location, color_location, texture_location;

    glfwSetErrorCallback(error_callback);

    if (!glfwInit())
        exit(EXIT_FAILURE);

    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 2);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 0);

    windows[0] = glfwCreateWindow(400, 400, "First", NULL, NULL);
    if (!windows[0])
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    glfwSetKeyCallback(windows[0], key_callback);

    glfwMakeContextCurrent(windows[0]);

    // Only enable vsync for the first of the windows to be swapped to
    // avoid waiting out the interval for each window
    glfwSwapInterval(1);

    // The contexts are created with the same APIs so the function
    // pointers should be re-usable between them
    gladLoadGLLoader((GLADloadproc)glfwGetProcAddress);

    // Create the OpenGL objects inside the first context, created above
    // All objects will be shared with the second context, created below
    {
        int x, y;
        char pixels[16 * 16];
        GLuint vertex_shader, fragment_shader;

        glGenTextures(1, &texture);
        glBindTexture(GL_TEXTURE_2D, texture);

        srand((unsigned int) glfwGetTimerValue());

        for (y = 0;  y < 16;  y++)
        {
            for (x = 0;  x < 16;  x++)
                pixels[y * 16 + x] = rand() % 256;
        }

        glTexImage2D(GL_TEXTURE_2D, 0, GL_RED, 16, 16, 0, GL_RED, GL_UNSIGNED_BYTE, pixels);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
        glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);

        vertex_shader = glCreateShader(GL_VERTEX_SHADER);
        glShaderSource(vertex_shader, 1, &vertex_shader_text, NULL);
        glCompileShader(vertex_shader);

        fragment_shader = glCreateShader(GL_FRAGMENT_SHADER);
        glShaderSource(fragment_shader, 1, &fragment_shader_text, NULL);
        glCompileShader(fragment_shader);

        program = glCreateProgram();
        glAttachShader(program, vertex_shader);
        glAttachShader(program, fragment_shader);
        glLinkProgram(program);

        mvp_location = glGetUniformLocation(program, "MVP");
        color_location = glGetUniformLocation(program, "color");
        texture_location = glGetUniformLocation(program, "texture");
        vpos_location = glGetAttribLocation(program, "vPos");

        glGenBuffers(1, &vertex_buffer);
        glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
        glBufferData(GL_ARRAY_BUFFER, sizeof(vertices), vertices, GL_STATIC_DRAW);
    }

    glUseProgram(program);
    glUniform1i(texture_location, 0);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);

    glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
    glEnableVertexAttribArray(vpos_location);
    glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE,
                          sizeof(vertices[0]), (void*) 0);

    windows[1] = glfwCreateWindow(400, 400, "Second", NULL, windows[0]);
    if (!windows[1])
    {
        glfwTerminate();
        exit(EXIT_FAILURE);
    }

    // Place the second window to the right of the first
    {
        int xpos, ypos, left, right, width;

        glfwGetWindowSize(windows[0], &width, NULL);
        glfwGetWindowFrameSize(windows[0], &left, NULL, &right, NULL);
        glfwGetWindowPos(windows[0], &xpos, &ypos);

        glfwSetWindowPos(windows[1], xpos + width + left + right, ypos);
    }

    glfwSetKeyCallback(windows[1], key_callback);

    glfwMakeContextCurrent(windows[1]);

    // While objects are shared, the global context state is not and will
    // need to be set up for each context

    glUseProgram(program);

    glEnable(GL_TEXTURE_2D);
    glBindTexture(GL_TEXTURE_2D, texture);

    glBindBuffer(GL_ARRAY_BUFFER, vertex_buffer);
    glEnableVertexAttribArray(vpos_location);
    glVertexAttribPointer(vpos_location, 2, GL_FLOAT, GL_FALSE,
                          sizeof(vertices[0]), (void*) 0);

    while (!glfwWindowShouldClose(windows[0]) &&
           !glfwWindowShouldClose(windows[1]))
    {
        int i;
        const vec3 colors[2] =
        {
            { 0.8f, 0.4f, 1.f },
            { 0.3f, 0.4f, 1.f }
        };

        for (i = 0;  i < 2;  i++)
        {
            int width, height;
            mat4x4 mvp;

            glfwGetFramebufferSize(windows[i], &width, &height);
            glfwMakeContextCurrent(windows[i]);

            glViewport(0, 0, width, height);

            mat4x4_ortho(mvp, 0.f, 1.f, 0.f, 1.f, 0.f, 1.f);
            glUniformMatrix4fv(mvp_location, 1, GL_FALSE, (const GLfloat*) mvp);
            glUniform3fv(color_location, 1, colors[i]);
            glDrawArrays(GL_TRIANGLE_FAN, 0, 4);

            glfwSwapBuffers(windows[i]);
        }

        glfwWaitEvents();
    }

    glfwTerminate();
    exit(EXIT_SUCCESS);
}
/*
// Add load texture function
int main(int argc, char **argv)
{
    sceneData = SceneData::GenerateDefaultScene();
    // setupCmdArgs(argc,argv);

    std::vector<GLFWwindow*> windows;

    const size_t RENDERER_COUNT = sceneData.cameraData.size();

    windows.resize(RENDERER_COUNT);

	// GLFWwindow* window;



	if (setupGLFW() < 0){
		return -1;
	}

    if (createGLFWwindow(windows[0], "Flocking") < 0){
        return -1;
    }

    glfwSetKeyCallback(windows[0], key_callback);
    // glfwSetCursorPosCallback(windows[0], mouse_callback);

    glfwMakeContextCurrent(windows[0]);

    // Only enable vsync for the first of the windows to be swapped to
    // avoid waiting out the interval for each window
    glfwSwapInterval(1);

    // The contexts are created with the same APIs so the function
    // pointers should be re-usable between them
	if (setupGLAD() < 0){
		return -1;
	}

    { // Initialize the scene here
        scene = Scene::GenerateDefaultScene(); // ground truth
        scene.RebindAllMeshes();
    }

    glEnable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);


    for (size_t windowIdx = 1; windowIdx < windows.size(); ++windowIdx){
        if (createGLFWwindow(windows[windowIdx], "Child", windows[0]) < 0){
            return -1;
        }

        // Place the second window to the right of the first
        {
            int xpos, ypos, left, right, width;

            glfwGetWindowSize(windows[windowIdx - 1], &width, NULL);
            glfwGetWindowFrameSize(windows[windowIdx - 1], &left, NULL, &right, NULL);
            glfwGetWindowPos(windows[windowIdx - 1], &xpos, &ypos);

            glfwSetWindowPos(windows[windowIdx], xpos + width + left + right, ypos);
        }

        glfwSetKeyCallback(windows[windowIdx], key_callback);
        glfwMakeContextCurrent(windows[windowIdx]);
        scene.RebindAllMeshes();
        glEnable(GL_DEPTH_TEST);
        glEnable(GL_BLEND);
        glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    }

    for (size_t rendererIdx = 0; rendererIdx < RENDERER_COUNT; ++rendererIdx){
        renderers.push_back(Renderer());

        // renderers[rendererIdx].mainCameraIdx = rendererIdx;
    }

    bool closingTime = false;

    while (!closingTime)
    {
        for (size_t windowIdx = 0; windowIdx < RENDERER_COUNT; ++windowIdx){
            int width, height;

            glfwGetFramebufferSize(windows[windowIdx], &width, &height);
            glfwMakeContextCurrent(windows[windowIdx]);

            glViewport(0, 0, width, height);

            // input
            // -----
            // processInput(windows[windowIdx], &renderer, &scene);
            renderers[windowIdx].Render(&scene);
            glfwSwapBuffers(windows[windowIdx]);
        }

        glfwWaitEvents();
        for (size_t windowIdx = 0; windowIdx < windows.size(); ++windowIdx){
            closingTime = closingTime || glfwWindowShouldClose(windows[windowIdx]);
        }

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

*/
// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------


void framebuffer_size_callback(GLFWwindow* window, int width, int height){
    //
}

float SPEED_MULTIPLIER = 10;
void processInput(GLFWwindow *window, Renderer* renderer, Scene* scene){ // abhorrent should only pass in renderer.
    // std::cout << "Frame Time: " << deltaTime <<  " ms" << std::endl;
    double mWidth = (double) SCR_WIDTH;
    double mHeight = (double) SCR_HEIGHT;

    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_NORMAL){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);
    }

    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }





    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
        scene->cameras[(renderer->mainCameraIdx + 1) % 2].moveForward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
        scene->cameras[(renderer->mainCameraIdx + 1) % 2].moveBackward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
        scene->cameras[(renderer->mainCameraIdx + 1) % 2].moveLeft(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
        scene->cameras[(renderer->mainCameraIdx + 1) % 2].moveRight(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS){
        scene->cameras[(renderer->mainCameraIdx + 1) % 2].moveUp(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
        scene->cameras[(renderer->mainCameraIdx + 1) % 2].moveDown(deltaTime);
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

