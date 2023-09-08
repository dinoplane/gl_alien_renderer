#include <glad/glad.h>
#include <GLFW/glfw3.h>

#include <glm/gtx/string_cast.hpp>
#include <shader_s.hpp>
#include <cube.hpp>
#include <pyramid.hpp>
#include <boid.hpp>


#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#include <iostream>

void framebuffer_size_callback(GLFWwindow* window, int width, int height);
void processInput(GLFWwindow *window);

// settings
const unsigned int SCR_WIDTH = 800;
const unsigned int SCR_HEIGHT = 600;


// camera
glm::vec3 cameraPos   = glm::vec3(0.0f, 0.0f,  0.5f);
glm::vec3 cameraFront = glm::vec3(0.0f, 0.0f, -1.0f);
glm::vec3 cameraUp    = glm::vec3(0.0f, 1.0f,  0.0f);

// timing
float deltaTime = 0.0f;	// time between current frame and last frame
float lastFrame = 0.0f;

bool firstMouse = true;

float yaw   = -90.0f;	// yaw is initialized to -90.0 degrees since a yaw of 0.0 results in a direction vector pointing to the right so we initially rotate a bit to the left.
float pitch =  0.0f;
float lastX =  800.0f / 2.0;
float lastY =  600.0 / 2.0;

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


	// // Compile shaders
	// unsigned int vertexShader = glCreateShader(GL_VERTEX_SHADER);
    // glShaderSource(vertexShader, 1, &vertexShaderSource, NULL);
    // glCompileShader(vertexShader);
    // // check for shader compile errors
    // int success;
    // char infoLog[512];
    // glGetShaderiv(vertexShader, GL_COMPILE_STATUS, &success);
    // if (!success)
    // {
    //     glGetShaderInfoLog(vertexShader, 512, NULL, infoLog);
    //     std::cout << "ERROR::SHADER::VERTEX::COMPILATION_FAILED\n" << infoLog << std::endl;
    // }

    // unsigned int fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
    // glShaderSource(fragmentShader, 1, &fragmentShaderSource, NULL);
    // glCompileShader(fragmentShader);
    // // check for shader compile errors
    // glGetShaderiv(fragmentShader, GL_COMPILE_STATUS, &success);
    // if (!success)
    // {
    //     glGetShaderInfoLog(fragmentShader, 512, NULL, infoLog);
    //     std::cout << "ERROR::SHADER::FRAGMENT::COMPILATION_FAILED\n" << infoLog << std::endl;
    // }


	// // Link shaders
	// shaderProgram = glCreateProgram();
    // glAttachShader(shaderProgram, vertexShader);
    // glAttachShader(shaderProgram, fragmentShader);
    // glLinkProgram(shaderProgram);
    // // check for linking errors
    // glGetProgramiv(shaderProgram, GL_LINK_STATUS, &success);
    // if (!success) {
    //     glGetProgramInfoLog(shaderProgram, 512, NULL, infoLog);
    //     std::cout << "ERROR::SHADER::PROGRAM::LINKING_FAILED\n" << infoLog << std::endl;
    // }
    // glDeleteShader(vertexShader);
    // glDeleteShader(fragmentShader);
}


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

// Create a camera class

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

    float sensitivity = 0.1f;
    xoffset *= sensitivity;
    yoffset *= sensitivity;

    yaw   += xoffset;
    pitch += yoffset;

    if(pitch > 89.0f)
        pitch = 89.0f;
    if(pitch < -89.0f)
        pitch = -89.0f;

    glm::vec3 direction;
    direction.x = cos(glm::radians(yaw)) * cos(glm::radians(pitch));
    direction.y = sin(glm::radians(pitch));
    direction.z = sin(glm::radians(yaw)) * cos(glm::radians(pitch));
    cameraFront = glm::normalize(direction);
}

// Add load texture function
int main()
{
	GLFWwindow* window;
	if (setupGLFW(window) < 0){
		return -1;
	}

	if (setupGLAD() < 0){
		return -1;
	}


    // glEnable(GL_DEPTH_TEST);

    // build and compile our shader program
    // ------------------------------------
    // unsigned int shaderProgram;
	// setupShaders(shaderProgram);
    Shader ourShader("./flocking/shader/exv.vs", "./flocking/shader/exf.fs");


    // uncomment this call to draw in wireframe polygons.
//    glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

    std::vector<unsigned int> textureID = loadAllTextures("sky.png", "awesomeface.png");
    // load and create a texture
    // -------------------------
    // unsigned int texture1, texture2;
    // // texture 1
    // // ---------
    // glGenTextures(1, &texture1);
    // glBindTexture(GL_TEXTURE_2D, texture1);
    //  // set the texture wrapping parameters
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);	// set texture wrapping to GL_REPEAT (default wrapping method)
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // // set texture filtering parameters
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // // load image, create texture and generate mipmaps
    // int width, height, nrChannels;
    // stbi_set_flip_vertically_on_load(true); // tell stb_image.h to flip loaded texture's on the y-axis.
    // // The FileSystem::getPath(...) is part of the GitHub repository so we can find files on any IDE/platform; replace it with your own image path.
    // unsigned char *data = stbi_load("sky.png", &width, &height, &nrChannels, 0);
    // if (data)
    // {
    //     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    //     glGenerateMipmap(GL_TEXTURE_2D);
    // }
    // else
    // {
    //     std::cout << "Failed to load texture" << std::endl;
    // }
    // stbi_image_free(data);
    // textureID.push_back(texture1);

    // // texture 2
    // // ---------
    // glGenTextures(1, &texture2);
    // glBindTexture(GL_TEXTURE_2D, texture2);
    // // set the texture wrapping parameters
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    // // set texture filtering parameters
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    // glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    // // load image, create texture and generate mipmaps
    // data = stbi_load("awesomeface.png", &width, &height, &nrChannels, 0);
    // if (data)
    // {
    //     // note that the awesomeface.png has transparency and thus an alpha channel, so make sure to tell OpenGL the data type is of GL_RGBA
    //     glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, width, height, 0, GL_RGBA, GL_UNSIGNED_BYTE, data);
    //     glGenerateMipmap(GL_TEXTURE_2D);
    // }
    // else
    // {
    //     std::cout << "Failed to load texture" << std::endl;
    // }
    // stbi_image_free(data);
    // textureID.push_back(texture2);

    // tell opengl for each sampler to which texture unit it belongs to (only has to be done once)
    // -------------------------------------------------------------------------------------------
    ourShader.use();
    ourShader.setInt("texture1", 0);
    ourShader.setInt("texture2", 1);




    // pass projection matrix to shader (as projection matrix rarely changes there's no need to do this per frame)
    // -----------------------------------------------------------------------------------------------------------
    ourShader.use();
    glm::mat4 view          = glm::mat4(1.0f); // make sure to initialize matrix to identity matrix first
    glm::mat4 projection    = glm::mat4(1.0f);
    projection = glm::perspective(glm::radians(60.0f), (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 2000.0f);
    ourShader.setMat4("projection", projection);
    // // camera/view transformation


    Cube test;
    test.init();
    test.scale(glm::vec3(100.0));


    Cube xcoord;
    xcoord.init();
    xcoord.translate(glm::vec3(20.0, 0.0, 0.0));
    xcoord.scale(glm::vec3(10.0));


    Cube zcoord;
    zcoord.init();
    zcoord.translate(glm::vec3(0.0, 0.0, 10.0));
    zcoord.scale(glm::vec3(10.0));

    std::cout << glm::to_string(test.modelMat) << std::endl;


    Boid b(glm::vec3(0.0, 3.0, 0.0));

    b.updatePosition();

    std::cout << glm::to_string(b.model.modelMat) << std::endl;
    glEnable(GL_DEPTH_TEST);
    glfwSetCursorPosCallback(window, mouse_callback);
    // render loop
    // -----------
    while (!glfwWindowShouldClose(window))
    {


        // input
        // -----
        processInput(window);

        // render
        // ------
        glClearColor(0.2f, 0.3f, 0.3f, 1.0f);
        glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

        // bind textures on corresponding texture units
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_2D, textureID[0]);
        glActiveTexture(GL_TEXTURE1);
        glBindTexture(GL_TEXTURE_2D, textureID[1]);

        // draw our first triangle

        view = glm::lookAt(cameraPos, cameraPos + cameraFront, cameraUp);

        ourShader.use();
        ourShader.setMat4("view", view);
        ourShader.setMat4("model", test.modelMat);
        test.render();

        ourShader.use();
        ourShader.setMat4("view", view);
        ourShader.setMat4("model", xcoord.modelMat);
        xcoord.render();

        ourShader.use();
        ourShader.setMat4("view", view);
        ourShader.setMat4("model", zcoord.modelMat);
        zcoord.render();



        ourShader.use();
        ourShader.setMat4("view", view);
        ourShader.setMat4("model", b.model.modelMat);
        b.model.render();

        b.updatePosition();

        // glBindVertexArray(VAO); // seeing as we only have a single VAO there's no need to bind it every time, but we'll do so to keep things a bit more organized
        // glDrawArrays(GL_TRIANGLES, 0, 3);
        // glDrawElements(GL_TRIANGLES, 6, GL_UNSIGNED_INT, 0);
        // glBindVertexArray(0); // no need to unbind it every time

        // glfw: swap buffers and poll IO events (keys pressed/released, mouse moved etc.)
        // -------------------------------------------------------------------------------
        glfwSwapBuffers(window);
        glfwPollEvents();
        deltaTime = glfwGetTime() - lastFrame;
        lastFrame = glfwGetTime();
    }

    // optional: de-allocate all resources once they've outlived their purpose:
    // ------------------------------------------------------------------------
    // glDeleteVertexArrays(1, &VAO);
    // glDeleteBuffers(1, &VBO);
    // glDeleteBuffers(1, &EBO);
    test.deleteBuffers();
    ourShader.deleteProgram();

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

    // GLFW_CURSOR_HIDDEN
    // glfwSetCursor
    // glfwdisa
    // glfwGetCursorPos(window, &mWidth, &mHeight);

    // if (glfwSetCursorPosCallback)
    if (glfwGetKey(window, GLFW_KEY_1) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_NORMAL){
        // std::cout << "MOOOSED " << deltaTime << std::endl;
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

        // glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_DISABLED);

    }

    if (glfwGetKey(window, GLFW_KEY_0) == GLFW_PRESS && glfwGetInputMode(window, GLFW_CURSOR) == GLFW_CURSOR_DISABLED){
        glfwSetInputMode(window, GLFW_CURSOR, GLFW_CURSOR_NORMAL);
    }


    // std::cout << mWidth << " " << mHeight << std::endl;
    if (glfwGetKey(window, GLFW_KEY_ESCAPE) == GLFW_PRESS)
        glfwSetWindowShouldClose(window, true);

    if (glfwGetKey(window, GLFW_KEY_W) == GLFW_PRESS){
        // std::cout << "PRESSED " << deltaTime << std::endl;
        cameraPos += cameraFront * deltaTime * SPEED_MULTIPLIER;
    }
    if (glfwGetKey(window, GLFW_KEY_S) == GLFW_PRESS){
        cameraPos -= cameraFront * deltaTime * SPEED_MULTIPLIER;
    }

    if (glfwGetKey(window, GLFW_KEY_A) == GLFW_PRESS){

        // std::cout << "PRESSED " << deltaTime << std::endl;
        cameraPos -= glm::normalize(glm::cross(cameraFront, cameraUp)) * deltaTime * SPEED_MULTIPLIER;
    }
    if (glfwGetKey(window, GLFW_KEY_D) == GLFW_PRESS){
        cameraPos += glm::normalize(glm::cross(cameraFront, cameraUp)) * deltaTime * SPEED_MULTIPLIER;
    }


}

// glfw: whenever the window size changed (by OS or user resize) this callback function executes
// ---------------------------------------------------------------------------------------------
void framebuffer_size_callback(GLFWwindow* window, int width, int height)
{
    // make sure the viewport matches the new window dimensions; note that width and
    // height will be significantly larger than specified on retina displays.
    glViewport(0, 0, width, height);
}

