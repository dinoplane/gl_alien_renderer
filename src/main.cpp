#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtx/string_cast.hpp>
#include <glm/gtc/random.hpp>

#include <fastgltf/core.hpp>
#include <fastgltf/types.hpp>
#include <fastgltf/tools.hpp>


#include <shader_s.hpp>
#include <camera.hpp>
#include <mesh.hpp>

#include <entity.hpp>
#include <scenedata.hpp>
#include <scene.hpp>
#include <cloth_system.hpp>
#include <renderer.hpp>
#include <scene_loader.hpp>
#include <model_loader.hpp>


#include <cstdlib>
#include <thread>
#include <barrier>
#include <chrono>
#include <cassert>

#include "stb_image.h"


#include <util.h>

#include <iostream>
#include <fstream>

#include "imgui.h"
#include "imgui_impl_glfw.h"
#include "imgui_impl_opengl3.h"
// #include <Tracy.hpp>




void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window, Renderer* renderer, Scene* scene);

// settings
const unsigned int SCR_WIDTH = 1000;
const unsigned int SCR_HEIGHT = 600;

const size_t RENDERER_SIDE = 3;
const size_t RENDERER_COUNT = RENDERER_SIDE * RENDERER_SIDE;

// int CHUNK_SIZE = std::max((NUM_BOIDS + (NUM_THREADS - 1)) / NUM_THREADS, (unsigned int) 1);
// // std::barrier sync_point(NUM_THREADS);

// // timing
double currentFrame = 0.0;
double deltaTime = 0.0;	// time between current frame and last frame
double lastFrame = 0.0;
double totalTime = 0.0;
uint32_t frameCount = 0;

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

uint32_t currRendererIdx = 0;
std::string scenePath; 
// glm::vec3 lightPos = glm::vec3(0.0, 15.0, 5.0);

fastgltf::Asset asset;

std::string cameraSettingsPath = "camera_settings.txt";

void setupCmdArgs(int argc, char **argv){
    if (argc == 2){
        scenePath = argv[1];
    } else scenePath = "./scene/minecraft.scn";

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
    // glfwSwapInterval(1); // Enable vsync
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

    if (key == GLFW_KEY_H && action == GLFW_PRESS){
        // Renderer::doCull = ~Renderer::doCull;
        Renderer::ToggleCull(scene);
        fmt::print("Culling: {}\n", Renderer::doCull);
    }

    if (key == GLFW_KEY_F && action == GLFW_PRESS){
        // Renderer::doCull = ~Renderer::doCull;
        Renderer::ToggleDebug();
    }

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
	// std::cout << src_str << ", " << type_str << ", " << severity_str << ", " << id << ": " << message << '\n';
}





// screen quad VAO

// void SetUpScreenQuad(GLuint* quadVAO, GLuint *quadVBO){
//     glCreateVertexArrays(1, quadVAO);
//     glEnableVertexArrayAttrib(*quadVAO, 0);
//     glEnableVertexArrayAttrib(*quadVAO, 1);

//     // position attribute
//     glVertexArrayAttribFormat(*quadVAO, 0, 2, GL_FLOAT, GL_FALSE, 0);
//     glVertexArrayAttribBinding(*quadVAO, 0, 0);

//     // texcoord attribute
//     glVertexArrayAttribFormat(*quadVAO, 1, 2, GL_FLOAT, GL_FALSE, 2 * sizeof(GL_FLOAT));
//     glVertexArrayAttribBinding(*quadVAO, 1, 0);


//     glCreateBuffers(1, quadVBO);
//     glNamedBufferStorage(*quadVBO, array_size(quadVertices) * sizeof(float) , &quadVertices, GL_DYNAMIC_STORAGE_BIT);

// }


static void RenderToFrame (const Scene& scene){
    // ZoneScoped;
    // First Pass
    for (const auto& particleSystem : scene.particleSystems) {
        particleSystem->CalculateForces();
    }
    for (Renderer& renderer : renderers){
        renderer.Render(scene);
    }

    // Second Pass
    glBindFramebuffer(GL_FRAMEBUFFER, 0);
    glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
    glClear(GL_COLOR_BUFFER_BIT);
    int startPixel = 0;
    int startPixelX = 0;
    
    int startPixelY = 0;
    for (size_t i = 0; i < RENDERER_SIDE; ++i){
        startPixelX = 0;
        for (size_t j = 0; j < RENDERER_SIDE; ++j){
            glBlitNamedFramebuffer(renderers[i * RENDERER_SIDE + j].FBO, 0,
             0, 0,
             renderers[i * RENDERER_SIDE + j].width, renderers[i * RENDERER_SIDE + j].height, 

            startPixelX, startPixelY,
            startPixelX + renderers[i * RENDERER_SIDE + j].width, startPixelY + renderers[i * RENDERER_SIDE + j].height,
            GL_COLOR_BUFFER_BIT, GL_LINEAR);
            // startPixel += renderers[i * RENDERER_SIDE + j].width;

            startPixelX += renderers[i * RENDERER_SIDE + j].width;
            
        }
        startPixelY += renderers[0].height;
    }

    // for (Renderer& renderer : renderers){

    //     // cout << "Rendering to screen " << startPixel << " " << renderer.width  << " " << renderer.height << endl;
    //     glBlitNamedFramebuffer(renderer.FBO, 0, 0, 0, renderer.width, renderer.height, startPixel, 0, startPixel + renderer.width, renderer.height, GL_COLOR_BUFFER_BIT, GL_LINEAR);
    //     startPixel += renderer.width;
    // }
// glPolygonMode( GL_FRONT_AND_BACK, GL_LINE );
    // screenShader->use();
    // glBindVertexArray(quadVAO);
    // glVertexArrayVertexBuffer(quadVAO, 0, quadVBO, 0, sizeof(float) * 4);
    // glDisable(GL_DEPTH_TEST);
    // glBindTexture(GL_TEXTURE_2D, renderer->FBOTexture);
    // glDrawArrays(GL_TRIANGLES, 0, 6);
// glPolygonMode( GL_FRONT_AND_BACK, GL_FILL );
}


/*
TODO List
- Add a way to load multiple textures for the indirect call
    - Procedurally create a texture atlas
    - Add appropriate material indices to the shader and model data
    - In the shader determine the formula for sampling the texture



*/



// Add load texture function
int main(int argc, char **argv)
{


	GLFWwindow* window;
	if (setupGLFW() < 0){
		return -1;
	}

    if (createGLFWwindow(window, "AlienGLRenderer") < 0){
        return -1;
    }

	if (setupGLAD() < 0){
		return -1;
	}

    setupCmdArgs(argc, argv);

    glfwSetCursorPosCallback(window, mouse_callback);
    glfwSetKeyCallback(window, key_callback);

    glfwMakeContextCurrent(window);

    glEnable(GL_DEBUG_OUTPUT);
    glEnable(GL_DEBUG_OUTPUT_SYNCHRONOUS);
    glDebugMessageCallback(message_callback, nullptr);
    glDebugMessageControl(GL_DONT_CARE, GL_DONT_CARE, GL_DEBUG_SEVERITY_NOTIFICATION, 0, nullptr, GL_FALSE);

    // Setup Dear ImGui context
    IMGUI_CHECKVERSION();
    ImGui::CreateContext();
    ImGuiIO& io = ImGui::GetIO(); (void)io;
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableKeyboard;     // Enable Keyboard Controls
    io.ConfigFlags |= ImGuiConfigFlags_NavEnableGamepad;      // Enable Gamepad Controls

    // Setup Dear ImGui style
    ImGui::StyleColorsDark();
    //ImGui::StyleColorsLight();

    // Setup Platform/Renderer backends
    ImGui_ImplGlfw_InitForOpenGL(window, true);
    const char* glsl_version = "#version 460";
    ImGui_ImplOpenGL3_Init(glsl_version);


    bool show_demo_window = true;
    bool show_another_window = false;
    ImVec4 clear_color = ImVec4(0.45f, 0.55f, 0.60f, 1.00f);

    currentFrame = glfwGetTime();
    deltaTime = currentFrame - lastFrame;
    lastFrame = currentFrame;

    // Setup Scene
    // fastgltf::Asset asset;
    // if (!ModelLoader::LoadGLTF("./resources/assets/models/teapot.gltf", &asset))
    // {
    //     assert(false);
    // }

    // for (auto& mesh : asset.meshes) {
    //     Primitive* teapot;

    //     if (!ModelLoader::LoadMesh(asset, mesh, teapot)) {
    //         assert(false);
    //     }
    // }
    //glfwSetWindowAttrib(window, GLFW_MAXIMIZED, GLFW_TRUE);



    SceneLoader::LoadScene(SceneLoader::LoadSceneData(scenePath), &scene);//Scene::GenerateDefaultScene();

    for ( const auto& [classname, entInstData] : scene.entityInstanceMap ){
        fmt::print("Material Index alpha cut off {}\n", entInstData.instModel.materials[0].alphaCutoff);
    }
    
    SceneLoader::LoadCameraSettings(cameraSettingsPath, &Renderer::allCameras, &Renderer::allDebugMeshes);

    for (size_t i = 0; i < RENDERER_COUNT; ++i){
        renderers.push_back(Renderer(SCR_WIDTH / RENDERER_SIDE, SCR_HEIGHT / RENDERER_SIDE, 0));
        renderers[i].shaderIdx = i;
    }
    Renderer::allCameras[0].setPerspectiveSize(SCR_WIDTH / RENDERER_SIDE, SCR_HEIGHT / RENDERER_SIDE);

    glfwMaximizeWindow(window);

    int width, height;
    glfwGetWindowSize(window, &width, &height);

    framebuffer_size_callback(window, width, height);

    fmt::print("Number of renderers: {}\n", renderers.size());
    while (!glfwWindowShouldClose(window))
    {
        // ZoneScoped;

        processInput(window, &renderers[currRendererIdx], &scene);
        RenderToFrame( scene );
    {
    // Start the Dear ImGui frame
        ImGui_ImplOpenGL3_NewFrame();
        ImGui_ImplGlfw_NewFrame();
        ImGui::NewFrame();

        // 1. Show the big demo window (Most of the sample code is in ImGui::ShowDemoWindow()! You can browse its code to learn more about Dear ImGui!).
        if (show_demo_window)
            ImGui::ShowDemoWindow(&show_demo_window);

        // 2. Show a simple window that we create ourselves. We use a Begin/End pair to create a named window.
        {
            ImGui::Begin("Hello, world!");                          // Create a window called "Hello, world!" and append into it.
            ImGui::Text("Culling %s", (Renderer::doCull ? "Enabled" : "Disabled"));
            ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", deltaTime * 1000.0f, 1.0f / deltaTime);

            int currRenderer = 0;
            int currCam = 0;

            for (Camera& cam : Renderer::allCameras) {
                std::string sepTextLabel = "Camera " + std::to_string(currCam);
                ImGui::SeparatorText(sepTextLabel.c_str());
                ImGui::Text("Position: %.2f, %.2f, %.2f\n", cam.position.x, cam.position.y, cam.position.z);
                ImGui::Text("Pitch: %.2f, Yaw: %.2f\n", cam.pitch, cam.yaw);
                ++currCam;
            }

            ImGui::End();
        }

        // 3. Show another simple window.
        if (show_another_window)
        {
            ImGui::Begin("Another Window", &show_another_window);   // Pass a pointer to our bool variable (the window will have a closing button that will clear the bool when clicked)
            ImGui::Text("Hello from another window!");
            if (ImGui::Button("Close Me"))
                show_another_window = false;
            ImGui::End();
        }

        // Rendering
        ImGui::Render();
        int display_w, display_h;
        glfwGetFramebufferSize(window, &display_w, &display_h);
        // glViewport(0, 0, display_w, display_h);
        // glClearColor(clear_color.x * clear_color.w, clear_color.y * clear_color.w, clear_color.z * clear_color.w, clear_color.w);
        // glClear(GL_COLOR_BUFFER_BIT);
        ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
        currentFrame = glfwGetTime();
        deltaTime = currentFrame - lastFrame;
        lastFrame = currentFrame;
        totalTime += deltaTime;
        ++frameCount;
    }
        glfwSwapBuffers(window);
        glfwPollEvents();
    }

    // glfw: terminate, clearing all previously allocated GLFW resources.
    // ------------------------------------------------------------------

    //ClothSystem* cl = static_cast<ClothSystem*> (scene.particleSystems[0].get());
    //std::cout << "Application force ms/frame " << cl->totalForceTime / frameCount * 1000 << " (" << cl->totalForceTime / totalTime * 100.0 << " %)" << std::endl;
    //std::cout << "Application solve ms/frame " << cl->totalSolveTime / frameCount * 1000 << " (" << cl->totalSolveTime / totalTime * 100.0 << " %)" << std::endl;
    //

    //std::cout << "Application average ms/frame " << totalTime / frameCount * 1000 << " (" << frameCount / totalTime << " FPS)" << std::endl;


    // Create and open a text file
    //std::ofstream MyFile("camera_settings.txt");

    //// Write to the file
    //MyFile << "Files can be tricky, but it is fun enough!";

    //// Close the file
    //MyFile.close();
    SceneLoader::SaveCameraSettings(cameraSettingsPath, Renderer::allCameras);



    ImGui_ImplOpenGL3_Shutdown();
    ImGui_ImplGlfw_Shutdown();
    ImGui::DestroyContext();

    glfwDestroyWindow(window);

    glfwTerminate();

   return 0;
}

// process all input: query GLFW whether relevant keys are pressed/released this frame and react accordingly
// ---------------------------------------------------------------------------------------------------------


void framebuffer_size_callback(GLFWwindow* window, int w, int h){
    // make sure the viewport matches the new window dimensions; note that w and
    // h will be significantly larger than specified on retina displays.
    glViewport(0, 0, w, h);
    
    for (Renderer& renderer : renderers){
        renderer.Resize(w / RENDERER_SIDE, h / RENDERER_SIDE);
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
    /*
    fmt::print("CamPos: {} {} {}\n", 
    Renderer::allCameras[renderer->mainCameraIdx].position.x, 
    Renderer::allCameras[renderer->mainCameraIdx].position.y, 
    Renderer::allCameras[renderer->mainCameraIdx].position.z);
    */
}

/*
Camera 0 Width: 960 Height: 1001 Yaw: -277.2 Pitch: 34.3 Fov: 30 Near: 0.1 Far: 30 Position: -1.82 1.51 -0.16
Camera 1 Width: 960 Height: 1001 Yaw: 817.9 Pitch: 46.2 Fov: 30 Near: 0.1 Far: 100 Position: 1.66 3.05 -1.33

*/