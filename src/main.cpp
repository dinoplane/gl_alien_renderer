#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/random.hpp>

#include <shader_s.hpp>
#include <camera.hpp>
#include <mesh.hpp>

#include <scene.hpp>
#include <entity.hpp>
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
uint currRendererIdx = 1;

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
int createGLFWwindow(GLFWwindow* &window, const char* windowName){
    // glfw window creation
    // --------------------

    window = glfwCreateWindow(SCR_WIDTH, SCR_HEIGHT, windowName, NULL, NULL);
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

static void error_callback(int error, const char* description)
{
    fprintf(stderr, "Error: %s\n", description);
}

int setupGLFW(){
	// glfw: initialize and configure
    // ------------------------------

    glfwSetErrorCallback(error_callback);

    glfwInit();
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 5);
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

void setupScreenShaders(){

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

        Renderer::allCameras[currRendererIdx].processMouseMovement(xoffset, yoffset);
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

        currRendererIdx = (currRendererIdx + 1) % renderers.size();
        std::cout << "Switching to renderer " << currRendererIdx << std::endl;
    }

    if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS){
        glfwSetWindowShouldClose(window, true);
    }
}


void RenderToFrame (Renderer* renderer, Scene* scene, Shader* screenShader, GLuint quadVAO, GLuint quadVBO){

    // First Pass
    renderer->Render(scene);

    // Second Pass
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
// glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    screenShader->use();
    glBindVertexArray(quadVAO);
    glVertexArrayVertexBuffer(quadVAO, 0, quadVBO, 0, sizeof(float) * 4);
    glDisable(GL_DEPTH_TEST);
    glBindTexture(GL_TEXTURE_2D, renderer->FBOTexture);
    glDrawArrays(GL_TRIANGLES, 0, 6);
// glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}

constexpr float quadVertices[] = { // vertex attributes for a quad that fills the entire screen in Normalized Device Coordinates.
    // positions   // texCoords
    -1.0f,  1.0f,  0.0f, 1.0f,
    -1.0f, -1.0f,  0.0f, 0.0f,
     1.0f, -1.0f,  1.0f, 0.0f,

    -1.0f,  1.0f,  0.0f, 1.0f,
     1.0f, -1.0f,  1.0f, 0.0f,
     1.0f,  1.0f,  1.0f, 1.0f
};

// screen quad VAO

void SetUpScreenQuad(GLuint* quadVAO, GLuint *quadVBO){
    glCreateVertexArrays(1, quadVAO);
    glEnableVertexArrayAttrib(*quadVAO, 0);
    glEnableVertexArrayAttrib(*quadVAO, 1);

    // position attribute
    glVertexArrayAttribFormat(*quadVAO, 0, 2, GL_FLOAT, GL_FALSE, 0);
    glVertexArrayAttribBinding(*quadVAO, 0, 0);

    // texcoord attribute
    glVertexArrayAttribFormat(*quadVAO, 1, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GL_FLOAT));
    glVertexArrayAttribBinding(*quadVAO, 1, 0);


    glCreateBuffers(1, quadVBO);
    glNamedBufferStorage(*quadVBO, array_size(quadVertices) * sizeof(float) , &quadVertices, GL_DYNAMIC_STORAGE_BIT);

}

void message_callback(GLenum source, GLenum type, GLuint id, GLenum severity, GLsizei length, GLchar const* message, void const* user_param)
{
	auto const src_str = [source]() {
		switch (source)
		{
		case GL_DEBUG_SOURCE_API: return "API";
		case GL_DEBUG_SOURCE_WINDOW_SYSTEM: return "WINDOW SYSTEM";
		case GL_DEBUG_SOURCE_SHADER_COMPILER: return "SHADER COMPILER";
		case GL_DEBUG_SOURCE_THIRD_PARTY: return "THIRD PARTY";
		case GL_DEBUG_SOURCE_APPLICATION: return "APPLICATION";
		case GL_DEBUG_SOURCE_OTHER: return "OTHER";
		}
	}();

	auto const type_str = [type]() {
		switch (type)
		{
		case GL_DEBUG_TYPE_ERROR: return "ERROR";
		case GL_DEBUG_TYPE_DEPRECATED_BEHAVIOR: return "DEPRECATED_BEHAVIOR";
		case GL_DEBUG_TYPE_UNDEFINED_BEHAVIOR: return "UNDEFINED_BEHAVIOR";
		case GL_DEBUG_TYPE_PORTABILITY: return "PORTABILITY";
		case GL_DEBUG_TYPE_PERFORMANCE: return "PERFORMANCE";
		case GL_DEBUG_TYPE_MARKER: return "MARKER";
		case GL_DEBUG_TYPE_OTHER: return "OTHER";
		}
	}();

	auto const severity_str = [severity]() {
		switch (severity) {
		case GL_DEBUG_SEVERITY_NOTIFICATION: return "NOTIFICATION";
		case GL_DEBUG_SEVERITY_LOW: return "LOW";
		case GL_DEBUG_SEVERITY_MEDIUM: return "MEDIUM";
		case GL_DEBUG_SEVERITY_HIGH: return "HIGH";
		}
	}();
	std::cout << src_str << ", " << type_str << ", " << severity_str << ", " << id << ": " << message << '\n';
}

// Add load texture function
int main(int argc, char **argv)
{
    const size_t RENDERER_COUNT = 2;

	GLFWwindow* window;
	if (setupGLFW() < 0){
		return -1;
	}

    if (createGLFWwindow(window, "Flocking") < 0){
        return -1;
    }

	if (setupGLAD() < 0){
		return -1;
	}

    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetKeyCallback(window, key_callback);

    glfwMakeContextCurrent(window);

glEnable(GL_DEBUG_OUTPUT);
glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
glDebugMessageCallback(message_callback, nullptr);
glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, nullptr, GL_FALSE);
    { // Initialize the scene here
        scene = Scene::GenerateDefaultScene();
    }

    for (size_t i = 0; i < RENDERER_COUNT; ++i){
        renderers.push_back(Renderer(800, 600));
    }
    Shader screenShader("./resources/shader/screen.vert", "./resources/shader/screen.frag");
    screenShader.use();

    uint quadVAO, quadVBO;
    SetUpScreenQuad(&quadVAO, &quadVBO);

    while (!glfwWindowShouldClose(window))
    {

        processInput(window, &renderers[currRendererIdx], &scene);
        // for (Renderer& renderer : renderers){
        //     renderer.Render(&scene);
        // }
        // renderers[currRendererIdx].Render(&scene);
        RenderToFrame(&renderers[currRendererIdx], &scene, &screenShader, quadVAO, quadVBO);


        glfwSwapBuffers(window);

        glfwPollEvents();

        deltaTime = glfwGetTime() - lastFrame;
        lastFrame = glfwGetTime();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------
    glfwTerminate();

   return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------


void framebuffer_size_callback(GLFWwindow* window, int width, int height){
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);

    for (Camera& camera : Renderer::allCameras){
        camera.setPerspectiveSize(width, height);
    }
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
        Renderer::allCameras[renderer->mainCameraIdx].moveForward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
        Renderer::allCameras[renderer->mainCameraIdx].moveBackward(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){
        Renderer::allCameras[renderer->mainCameraIdx].moveLeft(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
        Renderer::allCameras[renderer->mainCameraIdx].moveRight(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_SPACE) == GLFW_PRESS){
        Renderer::allCameras[renderer->mainCameraIdx].moveUp(deltaTime);
    }
    if (glfwGetKey(window, GLFW_KEY_LEFT_SHIFT) == GLFW_PRESS){
        Renderer::allCameras[renderer->mainCameraIdx].moveDown(deltaTime);
    }


}

