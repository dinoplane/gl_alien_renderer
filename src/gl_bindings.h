#ifndef GL_BINDINGS_H
#define GL_BINDINGS_H

#include <glad/glad.h>


// VERTEX SHADER INPUTS
const GLuint POSITION_ATTRIB_LOC = 0;
const GLuint NORMAL_ATTRIB_LOC = 1;
const GLuint TEXCOORD_ATTRIB_LOC = 2;
const GLuint MODEL_MATRIX_ATTRIB_LOC = 3;


// BASE_INSTANCE SHADERS BUFFERS
const GLuint PROJ_VIEW_UBO_BINDING = 0;
const GLuint MODEL_FROM_MESH_UBO_BINDING = 1;
const GLuint WORLD_FROM_MODEL_SSBO_BINDING = 2;
const GLuint INST_IS_RENDERED_SSBO_BINDING = 3;

// FRAGMENT SHADER INPUTS
const GLuint ALBEDO_TEXTURE_BINDING = 4;
const GLuint MATERIAL_UBO_BINDING = 5;

// CULLING SHADER BUFFERS
const GLuint FRUSTUM_CULL_DATA_UBO_BINDING = 4;


// SCREEN SHADER INPUTS
const GLuint SCREEN_TEXTURE_BINDING = 4;


#endif