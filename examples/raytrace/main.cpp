//
//   Copyright 2014 Pixar
//
//   Licensed under the Apache License, Version 2.0 (the "Apache License")
//   with the following modification; you may not use this file except in
//   compliance with the Apache License and the following modification to it:
//   Section 6. Trademarks. is deleted and replaced with:
//
//   6. Trademarks. This License does not grant permission to use the trade
//      names, trademarks, service marks, or product names of the Licensor
//      and its affiliates, except as required to comply with Section 4(c) of
//      the License and to reproduce the content of the NOTICE file.
//
//   You may obtain a copy of the Apache License at
//
//       http://www.apache.org/licenses/LICENSE-2.0
//
//   Unless required by applicable law or agreed to in writing, software
//   distributed under the Apache License with the above modification is
//   distributed on an "AS IS" BASIS, WITHOUT WARRANTIES OR CONDITIONS OF ANY
//   KIND, either express or implied. See the Apache License for the specific
//   language governing permissions and limitations under the Apache License.
//

#if defined(__APPLE__)
    #if defined(OSD_USES_GLEW)
        #include <GL/glew.h>
    #else
        #include <OpenGL/gl3.h>
    #endif
    #define GLFW_INCLUDE_GL3
    #define GLFW_NO_GLU
#else
    #include <stdlib.h>
    #include <GL/glew.h>
    #if defined(WIN32)
        #include <GL/wglew.h>
    #endif
#endif

#if defined(GLFW_VERSION_3)
    #include <GLFW/glfw3.h>
    GLFWwindow* g_window=0;
    GLFWmonitor* g_primary=0;
#else
    #include <GL/glfw.h>
#endif

#include <far/mesh.h>
#include <far/meshFactory.h>

#include <osd/error.h>
#include <osd/vertex.h>
#include <osd/glDrawContext.h>
#include <osd/glDrawRegistry.h>
#include <osd/glMesh.h>

#include <osd/cpuVertexBuffer.h>
#include <osd/cpuComputeContext.h>
#include <osd/cpuComputeController.h>

OpenSubdiv::OsdCpuVertexBuffer *g_cpuVertexBuffer = NULL;
OpenSubdiv::OsdCpuComputeContext *g_cpuComputeContext = NULL;

#include <common/shape_utils.h>
#include "../common/stopwatch.h"
#include "../common/simple_math.h"
#include "../common/gl_hud.h"
#include "../common/gl_common.h"

#include "scene.h"

#include <cfloat>
#include <vector>
#include <fstream>
#include <sstream>

static const char *s_VS =
    "#version 410\n"
    "in vec2 position;\n"
    "out vec2 uv;\n"
    "void main() {\n"
    "  gl_Position = vec4(position.x, position.y, 0, 1);\n"
    "  uv = (position+vec2(1))*0.5;\n"
    "}\n";
static const char *s_FS =
    "#version 410\n"
    "in vec2 uv;\n"
    "out vec4 outColor;\n"
    "uniform sampler2D tex;\n"
    "void main()\n"
    "{\n"
    "  vec2 texUV = textureSize(tex, 0)*uv;\n"
    "  vec4 c = vec4(0);\n"
    "  for(int i = 0; i < 8; ++i) {\n"
    "    for(int j = 0; j < i*2+1; ++j) {\n"
    "      ivec2 offset = min(ivec2(i*2-j, j), ivec2(i, i));\n"
    "      vec4 cs = texelFetch(tex, ivec2(texUV)-offset, 0);\n"
    "      c = mix(c, cs, 1-c.a);\n"
    "      c.a = max(c.a, cs.a);\n"
    "    }\n"
    "  }\n"
    "  outColor = c;\n"
    "}\n";

typedef OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>     OsdHbrMesh;
typedef OpenSubdiv::HbrVertex<OpenSubdiv::OsdVertex>   OsdHbrVertex;
typedef OpenSubdiv::HbrFace<OpenSubdiv::OsdVertex>     OsdHbrFace;
typedef OpenSubdiv::HbrHalfedge<OpenSubdiv::OsdVertex> OsdHbrHalfedge;

struct SimpleShape {
    std::string  name;
    Scheme       scheme;
    std::string  data;

    SimpleShape() { }
    SimpleShape( std::string const & idata, char const * iname, Scheme ischeme )
        : name(iname), scheme(ischeme), data(idata) { }
};

static void setCamera();

std::vector<SimpleShape> g_defaultShapes;
OpenSubdiv::FarMesh<OpenSubdiv::OsdVertex> *g_farMesh = NULL;

int g_currentShape = 0;

// GUI variables
int   g_displayStyle = Scene::SHADED,
      g_mbutton[3] = {0, 0, 0},
      g_partitioning = 1,
      g_running = 1;

float g_rotate[2] = {0, 0},
      g_dolly = 5,
      g_pan[2] = {0, 0},
      g_center[3] = {0, 0, 0},
      g_size = 0;

int   g_prev_x = 0,
      g_prev_y = 0;

int   g_frameBufferWidth = 1024,
    g_frameBufferHeight = 1024;

int   g_width = 1024,
      g_height = 1024;
std::vector<float> g_image;
int g_step = 8;
int g_stepIndex = 0;

GLhud g_hud;

// performance
float g_cpuTime = 0;
float g_gpuTime = 0;
Stopwatch g_fpsTimer;

// geometry
std::vector<float> g_orgPositions,
                   g_positions;

int g_level = 2;

struct Transform {
    float ModelViewMatrix[16];
    float ProjectionMatrix[16];
    float ModelViewProjectionMatrix[16];
} g_transformData;

GLuint g_vao = 0;
GLuint g_vbo = 0;
GLuint g_program = 0;

float g_eye[] = { 0, 0, 5, 1};
float g_lookat[] = {0, 0, 0, 1};
float g_up[] = {0, 1, 0, 0};

Scene g_scene;

//------------------------------------------------------------------------------
static void
initializeShapes( ) {

#include <shapes/catmark_cube_corner0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_corner0, "catmark_cube_corner0", kCatmark));

#include <shapes/catmark_cube_corner1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_corner1, "catmark_cube_corner1", kCatmark));

#include <shapes/catmark_cube_corner2.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_corner2, "catmark_cube_corner2", kCatmark));

#include <shapes/catmark_cube_corner3.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_corner3, "catmark_cube_corner3", kCatmark));

#include <shapes/catmark_cube_corner4.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_corner4, "catmark_cube_corner4", kCatmark));

#include <shapes/catmark_cube_creases0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_creases0, "catmark_cube_creases0", kCatmark));

#include <shapes/catmark_cube_creases1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube_creases1, "catmark_cube_creases1", kCatmark));

#include <shapes/catmark_cube.h>
    g_defaultShapes.push_back(SimpleShape(catmark_cube, "catmark_cube", kCatmark));

#include <shapes/catmark_dart_edgecorner.h>
    g_defaultShapes.push_back(SimpleShape(catmark_dart_edgecorner, "catmark_dart_edgecorner", kCatmark));

#include <shapes/catmark_dart_edgeonly.h>
    g_defaultShapes.push_back(SimpleShape(catmark_dart_edgeonly, "catmark_dart_edgeonly", kCatmark));

#include <shapes/catmark_edgecorner.h>
    g_defaultShapes.push_back(SimpleShape(catmark_edgecorner ,"catmark_edgecorner", kCatmark));

#include <shapes/catmark_edgeonly.h>
    g_defaultShapes.push_back(SimpleShape(catmark_edgeonly, "catmark_edgeonly", kCatmark));

#include <shapes/catmark_chaikin0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_chaikin0, "catmark_chaikin0", kCatmark));

#include <shapes/catmark_chaikin1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_chaikin1, "catmark_chaikin1", kCatmark));

#include <shapes/catmark_fan.h>
    g_defaultShapes.push_back(SimpleShape(catmark_fan, "catmark_fan", kCatmark));

#include <shapes/catmark_gregory_test1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_gregory_test1, "catmark_gregory_test1", kCatmark));

#include <shapes/catmark_gregory_test2.h>
    g_defaultShapes.push_back(SimpleShape(catmark_gregory_test2, "catmark_gregory_test2", kCatmark));

#include <shapes/catmark_gregory_test3.h>
    g_defaultShapes.push_back(SimpleShape(catmark_gregory_test3, "catmark_gregory_test3", kCatmark));

#include <shapes/catmark_gregory_test4.h>
    g_defaultShapes.push_back(SimpleShape(catmark_gregory_test4, "catmark_gregory_test4", kCatmark));

#include <shapes/catmark_hole_test1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_hole_test1, "catmark_hole_test1", kCatmark));

#include <shapes/catmark_hole_test2.h>
    g_defaultShapes.push_back(SimpleShape(catmark_hole_test2, "catmark_hole_test2", kCatmark));

#include <shapes/catmark_pyramid_creases0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_pyramid_creases0, "catmark_pyramid_creases0", kCatmark));

#include <shapes/catmark_pyramid_creases1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_pyramid_creases1, "catmark_pyramid_creases1", kCatmark));

#include <shapes/catmark_pyramid.h>
    g_defaultShapes.push_back(SimpleShape(catmark_pyramid, "catmark_pyramid", kCatmark));

#include <shapes/catmark_tent_creases0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_tent_creases0, "catmark_tent_creases0", kCatmark));

#include <shapes/catmark_tent_creases1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_tent_creases1, "catmark_tent_creases1", kCatmark));

#include <shapes/catmark_tent.h>
    g_defaultShapes.push_back(SimpleShape(catmark_tent, "catmark_tent", kCatmark));

#include <shapes/catmark_torus.h>
    g_defaultShapes.push_back(SimpleShape(catmark_torus, "catmark_torus", kCatmark));

#include <shapes/catmark_torus_creases0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_torus_creases0, "catmark_torus_creases0", kCatmark));

#include <shapes/catmark_square_hedit0.h>
    g_defaultShapes.push_back(SimpleShape(catmark_square_hedit0, "catmark_square_hedit0", kCatmark));

#include <shapes/catmark_square_hedit1.h>
    g_defaultShapes.push_back(SimpleShape(catmark_square_hedit1, "catmark_square_hedit1", kCatmark));

#include <shapes/catmark_square_hedit2.h>
    g_defaultShapes.push_back(SimpleShape(catmark_square_hedit2, "catmark_square_hedit2", kCatmark));

#include <shapes/catmark_square_hedit3.h>
    g_defaultShapes.push_back(SimpleShape(catmark_square_hedit3, "catmark_square_hedit3", kCatmark));

#include <shapes/catmark_square_hedit4.h>
    g_defaultShapes.push_back(SimpleShape(catmark_square_hedit4, "catmark_square_hedit4", kCatmark));

#include <shapes/catmark_bishop.h>
    g_defaultShapes.push_back(SimpleShape(catmark_bishop, "catmark_bishop", kCatmark));

#include <shapes/catmark_car.h>
    g_defaultShapes.push_back(SimpleShape(catmark_car, "catmark_car", kCatmark));

#include <shapes/catmark_helmet.h>
    g_defaultShapes.push_back(SimpleShape(catmark_helmet, "catmark_helmet", kCatmark));

#include <shapes/catmark_pawn.h>
    g_defaultShapes.push_back(SimpleShape(catmark_pawn, "catmark_pawn", kCatmark));

#include <shapes/catmark_rook.h>
    g_defaultShapes.push_back(SimpleShape(catmark_rook, "catmark_rook", kCatmark));
}

//------------------------------------------------------------------------------
static void
startRender() {
    g_image.clear();
    g_image.resize(g_width*g_height*4);
    g_stepIndex = g_step*g_step;
}

static void
setCamera() {
    // prepare view matrix
    double aspect = g_width/(double)g_height;
    identity(g_transformData.ModelViewMatrix);
    translate(g_transformData.ModelViewMatrix, g_pan[1], g_pan[0], -g_dolly);
    rotate(g_transformData.ModelViewMatrix, g_rotate[1], 0, 1, 0);
    rotate(g_transformData.ModelViewMatrix, g_rotate[0], 1, 0, 0);
    rotate(g_transformData.ModelViewMatrix, 90, 0, 0, 1);
    rotate(g_transformData.ModelViewMatrix, -90, 1, 0, 0);
    // translate(g_transformData.ModelViewMatrix,
    //           -g_center[0], -g_center[1], -g_center[2]);
    perspective(g_transformData.ProjectionMatrix,
                45.0f, (float)aspect, 0.01f, 500.0f);
    multMatrix(g_transformData.ModelViewProjectionMatrix,
               g_transformData.ModelViewMatrix,
               g_transformData.ProjectionMatrix);

    
    float invView[16];
    inverseMatrix(invView, g_transformData.ModelViewMatrix);

    g_eye[0] = g_eye[1] = g_eye[2] = 0;
    apply(g_eye, invView);

    g_lookat[0] = g_lookat[1] = 0;
    g_lookat[2] = -5;
    apply(g_lookat, invView);

    g_up[0] = g_up[2] = 1;
    g_up[1] = 0;
    apply(g_up, invView);


    startRender();
}

static void
updateGeom() {

    int nverts = (int)g_orgPositions.size() / 3;

    std::vector<float> vertex;
    vertex.reserve(nverts*3);

    const float *p = &g_orgPositions[0];

    for (int i = 0; i < nverts; ++i) {
        g_positions[i*3+0] = p[0];
        g_positions[i*3+1] = p[1];
        g_positions[i*3+2] = p[2];

        p += 3;
    }

    p = &g_orgPositions[0];
    const float *pp = &g_positions[0];
    for (int i = 0; i < nverts; ++i) {
        vertex.push_back(pp[0]);
        vertex.push_back(pp[1]);
        vertex.push_back(pp[2]);
        pp += 3;
    }

    g_cpuVertexBuffer->UpdateData(&vertex[0], 0, nverts);

    OpenSubdiv::OsdCpuComputeController controller;
    controller.Refine(g_cpuComputeContext,
                      g_farMesh->GetKernelBatches(),
                      g_cpuVertexBuffer);

    g_scene.Build(g_cpuVertexBuffer->BindCpuBuffer(),
                  g_cpuVertexBuffer->GetNumVertices(),
                  g_farMesh->GetPatchTables());
}

//------------------------------------------------------------------------------
static void
createOsdMesh( const std::string &shape, int level ){

    checkGLErrors("create osd enter");
    // generate Hbr representation from "obj" description
    OsdHbrMesh * hmesh = simpleHbr<OpenSubdiv::OsdVertex>(shape.c_str(), kCatmark, g_orgPositions);

    delete g_farMesh;
    delete g_cpuComputeContext;
    delete g_cpuVertexBuffer;

    // create farmesh
    OpenSubdiv::FarMeshFactory<OpenSubdiv::OsdVertex> meshFactory(hmesh, level, true);
    g_farMesh = meshFactory.Create();
    // Hbr mesh can be deleted
    delete hmesh;

    g_cpuComputeContext = OpenSubdiv::OsdCpuComputeContext::Create(
        g_farMesh->GetSubdivisionTables(), g_farMesh->GetVertexEditTables());
    g_cpuVertexBuffer = OpenSubdiv::OsdCpuVertexBuffer::Create(
        3, g_farMesh->GetNumVertices());

    // compute model bounding
    float min[3] = { FLT_MAX,  FLT_MAX,  FLT_MAX};
    float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
    for (size_t i=0; i <g_orgPositions.size()/3; ++i) {
        for(int j=0; j<3; ++j) {
            float v = g_orgPositions[i*3+j];
            min[j] = std::min(min[j], v);
            max[j] = std::max(max[j], v);
        }
    }
    for (int j=0; j<3; ++j) {
        g_center[j] = (min[j] + max[j]) * 0.5f;
        g_size += (max[j]-min[j])*(max[j]-min[j]);
    }
    g_size = sqrtf(g_size);

    g_positions.resize(g_orgPositions.size(),0.0f);

    updateGeom();
    setCamera();

    startRender();
}

//------------------------------------------------------------------------------
static void
fitFrame() {

    g_pan[0] = g_pan[1] = 0;
    g_dolly = g_size;
}

//------------------------------------------------------------------------------

static void
display() {

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glViewport(0, 0, g_frameBufferWidth, g_frameBufferHeight);

    static GLuint tex = 0;
    if (tex == 0) glGenTextures(1, &tex);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, tex);
    glPixelStorei(GL_UNPACK_ALIGNMENT, 1);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_NEAREST);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA,
                 g_width, g_height,
                 0, GL_RGBA, GL_FLOAT, &g_image[0]);

    glUseProgram(g_program);

    glBindVertexArray(g_vao);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE,
                          sizeof(GLfloat)*2, (void*)0);
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    glUseProgram(0);

    glBindTexture(GL_TEXTURE_2D, 0);
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    if (g_hud.IsVisible()) {
        g_fpsTimer.Stop();
        double fps = 1.0/g_fpsTimer.GetElapsed();
        g_fpsTimer.Start();

        if (g_stepIndex > 0) 
            g_hud.DrawString(10, -240, 1, 0, 0, "Rendering...");

        g_hud.DrawString(10, -180, "# of patches : %d", g_scene.GetNumPatches());
//        g_hud.DrawString(10, -60,  "GPU Draw   : %.3f ms", drawGpuTime);
//        g_hud.DrawString(10, -40,  "CPU Draw   : %.3f ms", drawCpuTime);
        g_hud.DrawString(10, -20,  "FPS        : %3.1f", fps);

        g_hud.Flush();
    }

    glFinish();

    glfwSwapBuffers(g_window);
    //checkGLErrors("display leave");
}

//------------------------------------------------------------------------------
static void
#if GLFW_VERSION_MAJOR>=3
motion(GLFWwindow *, double dx, double dy) {
    int x=(int)dx, y=(int)dy;
#else
motion(int x, int y) {
#endif

    if (g_mbutton[0] && !g_mbutton[1] && !g_mbutton[2]) {
        // orbit
        g_rotate[0] += x - g_prev_x;
        g_rotate[1] += y - g_prev_y;
        setCamera();
    } else if (!g_mbutton[0] && !g_mbutton[1] && g_mbutton[2]) {
        // pan
        g_pan[0] -= g_dolly*(x - g_prev_x)/g_width;
        g_pan[1] += g_dolly*(y - g_prev_y)/g_height;
        setCamera();
    } else if ((g_mbutton[0] && !g_mbutton[1] && g_mbutton[2]) or
               (!g_mbutton[0] && g_mbutton[1] && !g_mbutton[2])) {
        // dolly
        g_dolly -= g_dolly*0.01f*(x - g_prev_x);
        if(g_dolly <= 0.01) g_dolly = 0.01f;
        setCamera();
    }

    g_prev_x = x;
    g_prev_y = y;
}

//------------------------------------------------------------------------------
static void
#if GLFW_VERSION_MAJOR>=3
mouse(GLFWwindow *, int button, int state, int /* mods */) {
#else
mouse(int button, int state) {
#endif

    if (button == 0 && state == GLFW_PRESS && g_hud.MouseClick(g_prev_x, g_prev_y)) {
        display();
        return;
    }

    if (button < 3) {
        g_mbutton[button] = (state == GLFW_PRESS);
    }
}

//------------------------------------------------------------------------------
static void
uninitGL() {

    glDeleteVertexArrays(1, &g_vao);
    glDeleteBuffers(1, &g_vbo);

    delete g_farMesh;
    delete g_cpuVertexBuffer;
    delete g_cpuComputeContext;
}

//------------------------------------------------------------------------------
static void
#if GLFW_VERSION_MAJOR>=3
reshape(GLFWwindow *, int width, int height) {
#else
reshape(int width, int height) {
#endif

    g_width = width;
    g_height = height;

    startRender();

    g_frameBufferWidth = width;
    g_frameBufferHeight = height;

#if GLFW_VERSION_MAJOR>=3
    // window size might not match framebuffer size on a high DPI display
    glfwGetWindowSize(g_window, &g_width, &g_height);
#endif
    g_hud.Rebuild(g_width, g_height);
}

//------------------------------------------------------------------------------
#if GLFW_VERSION_MAJOR>=3
void windowClose(GLFWwindow*) {
    g_running = false;
}
#else
int windowClose() {
    g_running = false;
    return GL_TRUE;
}
#endif

//------------------------------------------------------------------------------
static void
#if GLFW_VERSION_MAJOR>=3
keyboard(GLFWwindow *, int key, int /* scancode */, int event, int /* mods */) {
#else
#define GLFW_KEY_ESCAPE GLFW_KEY_ESC
keyboard(int key, int event) {
#endif

    if (event == GLFW_RELEASE) return;
    if (g_hud.KeyDown(tolower(key))) return;

    switch (key) {
        case 'Q': g_running = 0; break;
        case 'F': fitFrame(); break;
        case GLFW_KEY_ESCAPE: g_hud.SetVisible(!g_hud.IsVisible()); break;
    }
}

//------------------------------------------------------------------------------
static void
rebuildOsdMesh()
{
    createOsdMesh( g_defaultShapes[ g_currentShape ].data, g_level );
}

static void
callbackLevel(int l)
{
    g_level = l;
    rebuildOsdMesh();
}

static void
callbackModel(int m)
{
    if (m < 0)
        m = 0;

    if (m >= (int)g_defaultShapes.size())
        m = (int)g_defaultShapes.size() - 1;

    g_currentShape = m;
    rebuildOsdMesh();
}

static void
callbackDisplayStyle(int b)
{
    g_displayStyle = b;
    g_scene.SetShadeMode((Scene::ShadeMode)b);

    startRender();
}

static void
initHUD()
{
#if GLFW_VERSION_MAJOR>=3
    // window size might not match framebuffer size on a high DPI display
    glfwGetWindowSize(g_window, &g_width, &g_height);
#endif
    g_hud.Init(g_width, g_height);

    int shading_pulldown = g_hud.AddPullDown("Shading (W)", 10, 10, 250, callbackDisplayStyle, 'w');
    g_hud.AddPullDownButton(shading_pulldown, "Shaded", Scene::SHADED, g_displayStyle==Scene::SHADED);
    g_hud.AddPullDownButton(shading_pulldown, "Ptex coord", Scene::PTEX_COORD, g_displayStyle==Scene::PTEX_COORD);

    for (int i = 1; i < 11; ++i) {
        char level[16];
        sprintf(level, "Lv. %d", i);
        g_hud.AddRadioButton(3, level, i==2, 10, 210+i*20, callbackLevel, i, '0'+(i%10));
    }

    int pulldown_handle = g_hud.AddPullDown("Shape (N)", -300, 10, 300, callbackModel, 'n');
    for (int i = 0; i < (int)g_defaultShapes.size(); ++i) {
        g_hud.AddPullDownButton(pulldown_handle, g_defaultShapes[i].name.c_str(),i);
    }   
}

//------------------------------------------------------------------------------
static void
initGL()
{
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    glGenVertexArrays(1, &g_vao);

    glGenBuffers(1, &g_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    static float pos[] = { -1, -1, 1, -1, -1,  1, 1,  1 };
    glBufferData(GL_ARRAY_BUFFER, sizeof(pos), pos, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, s_VS);
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, s_FS);
    g_program = glCreateProgram();
    glAttachShader(g_program, vertexShader);
    glAttachShader(g_program, fragmentShader);

    glLinkProgram(g_program);

    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
}

//------------------------------------------------------------------------------
static void
idle() {
    if (g_stepIndex == 0) return;

    double fov = 45.0f;

    //printf("%d x %d, %d\n", g_width, g_height, g_step);

    int s = 2*(g_step/2)*(g_step/2)+(g_step/2)+1;
    int index = (g_stepIndex*s)%(g_step*g_step);

    g_scene.Render(g_width, g_height, fov,
                   g_image,
                   g_eye, g_lookat, g_up, g_step, index);

    --g_stepIndex;
    display();
}

//------------------------------------------------------------------------------
static void
callbackError(OpenSubdiv::OsdErrorType err, const char *message)
{
    printf("OsdError: %d\n", err);
    printf("%s", message);
}

//------------------------------------------------------------------------------
static void
setGLCoreProfile()
{
#if GLFW_VERSION_MAJOR>=3
    #define glfwOpenWindowHint glfwWindowHint
    #define GLFW_OPENGL_VERSION_MAJOR GLFW_CONTEXT_VERSION_MAJOR
    #define GLFW_OPENGL_VERSION_MINOR GLFW_CONTEXT_VERSION_MINOR
#endif

    glfwOpenWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
#if not defined(__APPLE__)
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 4);
#ifdef OPENSUBDIV_HAS_GLSL_COMPUTE
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 3);
#else
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 2);
#endif
    
#else
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MAJOR, 3);
    glfwOpenWindowHint(GLFW_OPENGL_VERSION_MINOR, 2);
#endif
    glfwOpenWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
}

//------------------------------------------------------------------------------
int main(int argc, char ** argv)
{
    std::string str;
    for (int i = 1; i < argc; ++i) {
        if (!strcmp(argv[i], "-d"))
            g_level = atoi(argv[++i]);
        else {
            std::ifstream ifs(argv[1]);
            if (ifs) {
                std::stringstream ss;
                ss << ifs.rdbuf();
                ifs.close();
                str = ss.str();
                g_defaultShapes.push_back(SimpleShape(str.c_str(), argv[1], kCatmark));
            }
        }
    }
    initializeShapes();
    OsdSetErrorCallback(callbackError);

    if (not glfwInit()) {
        printf("Failed to initialize GLFW\n");
        return 1;
    }

    static const char windowTitle[] = "OpenSubdiv direct ray-tracing";

#define CORE_PROFILE
#ifdef CORE_PROFILE
    setGLCoreProfile();
#endif

#if GLFW_VERSION_MAJOR>=3
    if (not (g_window=glfwCreateWindow(g_width, g_height, windowTitle, NULL, NULL))) {
        printf("Failed to open window.\n");
        glfwTerminate();
        return 1;
    }
    glfwMakeContextCurrent(g_window);

    // accommocate high DPI displays (e.g. mac retina displays)
    glfwGetFramebufferSize(g_window, &g_frameBufferWidth, &g_frameBufferHeight);
    glfwSetFramebufferSizeCallback(g_window, reshape);

    glfwSetKeyCallback(g_window, keyboard);
    glfwSetCursorPosCallback(g_window, motion);
    glfwSetMouseButtonCallback(g_window, mouse);
    glfwSetWindowCloseCallback(g_window, windowClose);
#else
    if (glfwOpenWindow(g_width, g_height, 8, 8, 8, 8, 24, 8, GLFW_WINDOW) == GL_FALSE) {
        printf("Failed to open window.\n");
        glfwTerminate();
        return 1;
    }
    glfwSetWindowTitle(windowTitle);
    glfwSetKeyCallback(keyboard);
    glfwSetMousePosCallback(motion);
    glfwSetMouseButtonCallback(mouse);
    glfwSetWindowSizeCallback(reshape);
    glfwSetWindowCloseCallback(windowClose);
#endif


#if defined(OSD_USES_GLEW)
#ifdef CORE_PROFILE
    // this is the only way to initialize glew correctly under core profile context.
    glewExperimental = true;
#endif
    if (GLenum r = glewInit() != GLEW_OK) {
        printf("Failed to initialize glew. Error = %s\n", glewGetErrorString(r));
        exit(1);
    }
#ifdef CORE_PROFILE
    // clear GL errors which was generated during glewInit()
    glGetError();
#endif
#endif

    initGL();

    glfwSwapInterval(0);

    initHUD();
    rebuildOsdMesh();

    g_image.resize(g_width*g_height*4);

    while (g_running) {
        idle();

        glfwPollEvents();

        if (g_stepIndex == 0) glfwWaitEvents();
    }

    uninitGL();
    glfwTerminate();
}

//------------------------------------------------------------------------------
