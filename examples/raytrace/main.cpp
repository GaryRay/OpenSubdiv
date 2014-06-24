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

#define NEED_HBR_FACE_INDEX

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
    "  uv = (-position.yx+vec2(1))*0.5;\n"
    "}\n";
static const char *s_FS =
    "#version 410\n"
    "in vec2 uv;\n"
    "out vec4 outColor;\n"
    "uniform sampler2D tex;\n"
    "void main()\n"
    "{\n"
    "  ivec2 texUV = ivec2(textureSize(tex, 0)*uv);\n"
    "  vec4 c = vec4(0.1, 0.1, 0.1, 0);\n"
    "  vec4 cs = texelFetch(tex, texUV, 0);\n"
    "  vec4 c0 = texelFetch(tex, ivec2(texUV.x & ~7, texUV.y & ~7), 0);\n"
    "  vec4 c1 = texelFetch(tex, ivec2(texUV.x & ~3, texUV.y & ~3), 0);\n"
    "  vec4 c2 = texelFetch(tex, ivec2(texUV.x & ~1, texUV.y & ~1), 0);\n"
    "  c = mix(c, c0, c0.a);\n"
    "  c = mix(c, c1, c1.a);\n"
    "  c = mix(c, c2, c2.a);\n"
    "  c = mix(c, cs, cs.a);\n"
    "  outColor = c;\n"
    "}\n";
static const char *s0_FS =
    "#version 410\n"
    "in vec2 uv;\n"
    "out vec4 outColor;\n"
    "uniform sampler2D tex;\n"
    "void main()\n"
    "{\n"
    "  vec2 texUV = textureSize(tex, 0)*uv;\n"
    "  vec4 c = vec4(0.1, 0.1, 0.1, 0);\n"
    "  vec4 cs = texelFetch(tex, ivec2(texUV), 0);\n"
    "  outColor = mix(c, cs, cs.a);\n"
    "}\n";

static const char *s_VS_BVH =
    "#version 410\n"
    "in vec3 minPos;\n"
    "in vec3 maxPos;\n"
    "in int flag;\n"
    "out vec3 gMinPos;\n"
    "out vec3 gMaxPos;\n"
    "out vec4 gColor;\n"
    "void main() {\n"
    "  gMinPos = minPos;\n"
    "  gMaxPos = maxPos;\n"
    "  gColor = (flag==1) ? vec4(0.0, 0.8, 0.0, 0.5) : vec4(0.1, 0.3, 0.1, 0.5);\n"
    "}\n";

static const char *s_GS_BVH =
    "#version 410\n"
    "layout(points) in;\n"
    "layout(line_strip, max_vertices = 24) out;\n"
    "in vec3 gMinPos[];\n"
    "in vec3 gMaxPos[];\n"
    "in vec4 gColor[];\n"
    "out vec4 fColor;\n"
    "uniform mat4 ModelViewProjectionMatrix;\n"
    "void emit(vec3 p0, vec3 p1){\n"
    "  gl_Position = ModelViewProjectionMatrix * vec4(p0, 1);\n"
    "  EmitVertex();\n"
    "  gl_Position = ModelViewProjectionMatrix * vec4(p1, 1);\n"
    "  EmitVertex();\n"
    "  EndPrimitive();\n"
    "}\n"
    "void main() {\n"
    "  vec3 p0 = gMinPos[0];\n"
    "  vec3 p1 = gMaxPos[0];\n"
    "  fColor = gColor[0];\n"
    "  emit(vec3(p0.x, p0.y, p0.z), vec3(p1.x, p0.y, p0.z));\n"
    "  emit(vec3(p0.x, p0.y, p1.z), vec3(p1.x, p0.y, p1.z));\n"
    "  emit(vec3(p0.x, p1.y, p0.z), vec3(p1.x, p1.y, p0.z));\n"
    "  emit(vec3(p0.x, p1.y, p1.z), vec3(p1.x, p1.y, p1.z));\n"
    "  emit(vec3(p0.x, p0.y, p0.z), vec3(p0.x, p1.y, p0.z));\n"
    "  emit(vec3(p0.x, p0.y, p1.z), vec3(p0.x, p1.y, p1.z));\n"
    "  emit(vec3(p1.x, p0.y, p0.z), vec3(p1.x, p1.y, p0.z));\n"
    "  emit(vec3(p1.x, p0.y, p1.z), vec3(p1.x, p1.y, p1.z));\n"
    "  emit(vec3(p0.x, p0.y, p0.z), vec3(p0.x, p0.y, p1.z));\n"
    "  emit(vec3(p0.x, p1.y, p0.z), vec3(p0.x, p1.y, p1.z));\n"
    "  emit(vec3(p1.x, p0.y, p0.z), vec3(p1.x, p0.y, p1.z));\n"
    "  emit(vec3(p1.x, p1.y, p0.z), vec3(p1.x, p1.y, p1.z));\n"
    "}\n";

static const char *s_FS_BVH =
    "#version 410\n"
    "in vec4 fColor;\n"
    "out vec4 outColor;\n"
    "void main() {\n"
    "  outColor = fColor;\n"
    "}\n";

float g_debugScale = 0.01f;
static const char *s_VS_Debug =
    "#version 410\n"
    "in vec2 position;\n"
    "out vec2 uv;\n"
    "uniform float debugScale=0.01;\n"
    "void main() {\n"
    "  vec2 pos = position.xy * 0.25 + vec2(0.75, -0.75);\n"
    "  uv = (-position.yx*debugScale+vec2(1))*0.5;\n"
    "  gl_Position = vec4(pos.x, pos.y, 0, 1);\n"
    "}\n";

typedef OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>     OsdHbrMesh;
typedef OpenSubdiv::HbrVertex<OpenSubdiv::OsdVertex>   OsdHbrVertex;
typedef OpenSubdiv::HbrFace<OpenSubdiv::OsdVertex>     OsdHbrFace;
typedef OpenSubdiv::HbrHalfedge<OpenSubdiv::OsdVertex> OsdHbrHalfedge;

enum HudCheckBox { kHUD_CB_DISPLAY_BVH,
                   kHUD_CB_BLOCK_FILL,
                   kHUD_CB_WATERTIGHT,
                   kHUD_CB_CROP_UV,
                   kHUD_CB_BEZIER_CLIP,
                   kHUD_CB_PRE_TESSELLATE,
                   kHUD_CB_ANIMATE,
                   kHUD_CB_DEBUG };

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
OsdHbrMesh *g_hbrMesh = NULL;

int g_currentShape = 0;

// GUI variables
int   g_displayStyle = Scene::SHADED,
      g_drawBVH = false,
      g_blockFill = true,
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
Stopwatch g_renderTimer;
float g_hbrTime = 0;
float g_farTime = 0;
float g_subdivTime = 0;
float g_convertTime = 0;
float g_tessellateTime = 0;
float g_bvhTime = 0;
float g_renderTime = 0;

// geometry
std::vector<float> g_orgPositions;

int g_level = 2;
int g_preTess = 0;
int g_preTessLevel = 1;
int g_intersectKernel = 1;
int g_watertight = 1;
int g_cropUV = 1;
int g_bezierClip = 1;
int g_debug = 0;
float g_uvMargin = 0.01f;
float g_displaceScale = 0.0f;
float g_displaceFreq = 100.0f;

int g_animate = 0;
int g_frame = 0;

struct Transform {
    float ModelViewMatrix[16];
    float ProjectionMatrix[16];
    float ModelViewProjectionMatrix[16];
} g_transformData;

GLuint g_vao = 0;
GLuint g_vaoBVH = 0;
GLuint g_vbo = 0;
GLuint g_programSimpleFill = 0;
GLuint g_programBlockFill = 0;
GLuint g_programBVH = 0;
GLuint g_programDebug = 0;

float g_eye[] = { 0, 0, 5, 1};
float g_lookat[] = {0, 0, 0, 1};
float g_up[] = {0, 1, 0, 0};

Scene g_scene;
std::vector<int> g_vertexParentIDs;
std::vector<int> g_farToHbrVertexRemap;

//------------------------------------------------------------------------------
static void
initializeShapes( ) {

#include <shapes/catmark_extraordinary.h>
    g_defaultShapes.push_back(SimpleShape(catmark_extraordinary, "catmark_extraordinary", kCatmark));

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

    g_renderTimer.Start();
}

static void
setCamera() {
    // prepare view matrix
    double aspect = g_width/(double)g_height;
    identity(g_transformData.ModelViewMatrix);
#if 0
    translate(g_transformData.ModelViewMatrix, g_pan[1], g_pan[0], -g_dolly);
    rotate(g_transformData.ModelViewMatrix, g_rotate[1], 0, 1, 0);
    rotate(g_transformData.ModelViewMatrix, g_rotate[0], 1, 0, 0);
    rotate(g_transformData.ModelViewMatrix, 90, 0, 0, 1);
    rotate(g_transformData.ModelViewMatrix, -90, 1, 0, 0);
#else
    translate(g_transformData.ModelViewMatrix, -g_pan[0], -g_pan[1], -g_dolly);
    rotate(g_transformData.ModelViewMatrix, g_rotate[1], 1, 0, 0);
    rotate(g_transformData.ModelViewMatrix, g_rotate[0], 0, 1, 0);
    rotate(g_transformData.ModelViewMatrix, -90, 1, 0, 0);
#endif
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
dumpCamera()
{
  const char *filename = "camera.dat";

  FILE* fp = fopen(filename, "w");
  if (!fp) {
    fprintf(stderr, "Failed tp open file: %s\n", filename);
    return;
  }

  fprintf(fp, "%f %f %f\n", g_pan[0], g_pan[1], g_dolly);
  fprintf(fp, "%f %f\n", g_rotate[0], g_rotate[1]);

  fclose(fp);

  printf("camera data dumped.\n"); fflush(stdout);
}

static void
loadCamera()
{
  const char *filename = "camera.dat";
  FILE* fp = fopen(filename, "r");
  if (!fp) {
    fprintf(stderr, "Failed tp open file: %s\n", filename);
    return;
  }

  fscanf(fp, "%f %f %f\n", &g_pan[0], &g_pan[1], &g_dolly);
  fscanf(fp, "%f %f\n", &g_rotate[0], &g_rotate[1]);
  fclose(fp);
  printf("camera data loaded.\n"); fflush(stdout);
}

static void
updateGeom() {

    int nverts = (int)g_orgPositions.size() / 3;

    std::vector<float> vertex;
    vertex.reserve(nverts*3);

    float r = g_animate ? sin(g_frame*0.1f) : 0.0f;
    const float *pp = &g_orgPositions[0];
    for (int i = 0; i < nverts; ++i) {
        float x = pp[0];
        float y = pp[1];
        float z = pp[2];
        float ct = cos(y * r);
        float st = sin(y * r);
        vertex.push_back(x + ct + z * st);
        vertex.push_back(y);
        vertex.push_back(- x * st + z * ct);
        pp += 3;
    }

    g_cpuVertexBuffer->UpdateData(&vertex[0], 0, nverts);

    Stopwatch s;
    s.Start();
    OpenSubdiv::OsdCpuComputeController controller;
    controller.Refine(g_cpuComputeContext,
                      g_farMesh->GetKernelBatches(),
                      g_cpuVertexBuffer);
    s.Stop();
    g_subdivTime = s.GetElapsed() * 1000.0f;

    g_scene.SetWatertight(g_watertight);

    s.Start();
    g_scene.BezierConvert(g_cpuVertexBuffer->BindCpuBuffer(),
                          g_cpuVertexBuffer->GetNumVertices(),
                          g_farMesh->GetPatchTables(),
                          //g_vertexParentIDs,
                          g_farToHbrVertexRemap,
                          g_hbrMesh,
                          g_displaceScale/*bound*/);
    s.Stop();
    g_convertTime = s.GetElapsed() * 1000.0f;

    if (g_preTess) {
        s.Start();
        g_scene.Tessellate(g_preTessLevel);
        s.Stop();
        g_tessellateTime = s.GetElapsed() * 1000.0f;
    } else {
        g_tessellateTime = 0;
    }

    s.Start();
    g_scene.Build();
    s.Stop();

    g_bvhTime = s.GetElapsed() * 1000.0f;

    g_scene.VBOBuild();

    startRender();
}

//------------------------------------------------------------------------------
static void
createOsdMesh( const std::string &shape, int level ){

    checkGLErrors("create osd enter");
    // generate Hbr representation from "obj" description

    Stopwatch s;

    s.Start();
    if (g_hbrMesh) delete g_hbrMesh;
    g_hbrMesh = simpleHbr<OpenSubdiv::OsdVertex>(shape.c_str(), kCatmark, g_orgPositions);
    s.Stop();

    g_hbrTime = s.GetElapsed() * 1000.0f;

    delete g_farMesh;
    delete g_cpuComputeContext;
    delete g_cpuVertexBuffer;

    // create farmesh
    s.Start();
    OpenSubdiv::FarMeshFactory<OpenSubdiv::OsdVertex> meshFactory(g_hbrMesh, level, true);
    g_farMesh = meshFactory.Create();
    s.Stop();


    // create index reduction table.
    std::vector<int> remap = meshFactory.GetRemappingTable();
    // for (int i = 0; i < remap.size(); ++i) {
    //     printf("%d ", remap[i]);
    // }
//    printf("\n");

    std::vector<int> parentIDs(remap.size());
    std::vector<int> farToHbrVertexRemap(remap.size());
    for (int i = 0; i < (int)remap.size(); ++i) {
        OsdHbrVertex *vertex = g_hbrMesh->GetVertex(i);
        int parentID = i;
        do {
            parentID = vertex->GetID();
            vertex = vertex->GetParentVertex();
        } while(vertex);

        parentIDs[remap[i]] = remap[parentID];
        farToHbrVertexRemap[remap[i]] = i;
    }
    // for (int i = 0; i < parentIDs.size(); ++i) {
    //     printf("%d ", parentIDs[i]);
    // }
    // printf("\n");
    g_vertexParentIDs = parentIDs;
    g_farToHbrVertexRemap = farToHbrVertexRemap;


    g_farTime = s.GetElapsed() * 1000.0f;

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

    setCamera();
    updateGeom();
}

float g_selectedColor[3];
static void
debugTrace(int x, int y) {
    float *p = &g_image[(y*g_width + x)*4];
    g_selectedColor[0] = p[0];
    g_selectedColor[1] = p[1];
    g_selectedColor[2] = p[2];
}

//------------------------------------------------------------------------------
static void
fitFrame() {

    g_pan[0] = g_pan[1] = 0;
    g_dolly = g_size;

    setCamera();
}

//------------------------------------------------------------------------------
static void
displayBVH() {

    glDisable(GL_DEPTH_TEST);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBindVertexArray(g_vaoBVH);

    glBindBuffer(GL_ARRAY_BUFFER, g_scene.GetVBO());
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(BVHNode), (void*)0);
    glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(BVHNode), (void*)(sizeof(float)*3));
    glVertexAttribIPointer(2, 1, GL_INT, sizeof(BVHNode), (void*)(sizeof(float)*6));

    glUseProgram(g_programBVH);

    glUniformMatrix4fv(0, 1, GL_FALSE, g_transformData.ModelViewProjectionMatrix);

    glDrawArrays(GL_POINTS, 0, g_scene.GetNumBVHNode());

    glUseProgram(0);

    glBindVertexArray(0);

    glDisable(GL_BLEND);

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

    GLuint program = g_blockFill ? g_programBlockFill : g_programSimpleFill;
    glUseProgram(program);

    glBindVertexArray(g_vao);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    glEnableVertexAttribArray(0);
    glVertexAttribPointer(0, 2, GL_FLOAT, GL_FALSE,
                          sizeof(GLfloat)*2, (void*)0);

    //int loc = glGetUniformLocation(program, "scale");
    glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);

    if (g_debug) {
        glUseProgram(g_programDebug);
        glDrawArrays(GL_TRIANGLE_STRIP, 0, 4);
    }

    glUseProgram(0);

    glBindTexture(GL_TEXTURE_2D, 0);
    glDisableVertexAttribArray(0);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
    glBindVertexArray(0);

    if (g_drawBVH) {
        displayBVH();
    }

    if (g_hud.IsVisible()) {
        g_fpsTimer.Stop();
        double fps = 1.0/g_fpsTimer.GetElapsed();
        g_fpsTimer.Start();

        g_hud.DrawString(10, -300, "# of patches      : %8d", g_scene.GetNumPatches());
        g_hud.DrawString(10, -280, "PreTess lv        : %d (+/-)", g_preTessLevel);
        g_hud.DrawString(10, -260, "# of tris         : %8d", g_scene.GetNumTriangles());

        g_hud.DrawString(10, -220, "Memory            : %8.1f MB", g_scene.GetMemoryUsage()/1024.0/1024.0);

        g_hud.DrawString(10, -200, "Hbr time          : %8.1f ms", g_hbrTime);
        g_hud.DrawString(10, -180, "Far time          : %8.1f ms", g_farTime);
        g_hud.DrawString(10, -160, "Subdivision time  : %8.1f ms", g_subdivTime);

        g_hud.DrawString(10, -140, "Bspline to Bezier : %8.1f ms", g_convertTime);
        g_hud.DrawString(10, -120, "Bezier to Tri     : %8.1f ms", g_tessellateTime);
        g_hud.DrawString(10,  -80, "BVH build         : %8.1f ms", g_bvhTime);
        if (g_renderTime > 0) {
            g_hud.DrawString(10, -40,  "Render time       : %5.3f s", g_renderTime);
        } else {
            g_hud.DrawString(10, -40,  1, 0, 0, "Render time       : %2.0f%%",
                             100*(1-g_stepIndex/float(g_step*g_step)));
        }
        g_hud.DrawString(10, -20,  "FPS               : %3.1f", fps);


        if (g_debug) {
            g_hud.DrawString(-300, -300,
                             g_selectedColor[0], g_selectedColor[1], g_selectedColor[2],
                             "COLOR : %f %f %f", 
                             g_selectedColor[0], g_selectedColor[1], g_selectedColor[2]);
            g_hud.DrawString(g_width/2-10, g_height/2-5, "[ ]");
        }

        g_hud.Flush();
    }


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

    if (g_hud.MouseCapture()) {
        // check gui (for slider)
        g_hud.MouseMotion(x, y);
    } else if (g_debug) {
        //
    } else if (g_mbutton[0] && !g_mbutton[1] && !g_mbutton[2]) {
        // orbit
        g_rotate[0] += (x - g_prev_x);
        g_rotate[1] += (y - g_prev_y);
        setCamera();
    } else if (!g_mbutton[0] && !g_mbutton[1] && g_mbutton[2]) {
        // pan
        g_pan[0] -= 0.5*g_dolly*(x - g_prev_x)/g_width;
        g_pan[1] += 0.5*g_dolly*(y - g_prev_y)/g_height;
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
mouse(GLFWwindow *, int button, int state, int mods) {
#else
mouse(int button, int state) {
#endif

    if (state == GLFW_RELEASE) {
        g_hud.MouseRelease();
        g_mbutton[0] = g_mbutton[1] = g_mbutton[2] = 0;
    }

    if (button == 0 && state == GLFW_PRESS && g_hud.MouseClick(g_prev_x, g_prev_y)) {
        display();
        return;
    }

    if (mods & GLFW_MOD_SHIFT) {
        // for mac
        button = 2;
    }
    if (button < 3) {
        g_mbutton[button] = (state == GLFW_PRESS);
    }

    if (g_debug && button == 0) {
        float x = g_prev_x/float(g_width);
        float y = g_prev_y/float(g_height);
        x = (((x - 0.75)/0.25)*2 - 1);
        y = (((y - 0.75)/0.25)*2 - 1);
        if (x >= -1 and y >=-1 and x <= 1 and y <= 1) {
            x = g_width/2 - x*g_width*g_debugScale*0.5;
            y = g_height/2 + y*g_height*g_debugScale*0.5;

            debugTrace((int)y, (int)x);
            display();
        }
    }

}

//------------------------------------------------------------------------------
static void
uninitGL() {

    glDeleteVertexArrays(1, &g_vao);
    glDeleteVertexArrays(1, &g_vaoBVH);
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

    g_frameBufferWidth = width;
    g_frameBufferHeight = height;

#if GLFW_VERSION_MAJOR>=3
    // window size might not match framebuffer size on a high DPI display
    glfwGetWindowSize(g_window, &g_width, &g_height);
#endif
    g_hud.Rebuild(g_width, g_height);

    setCamera();
    startRender();
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

#if GLFW_VERSION_MAJOR<3
#error "please use glfw3."
#endif

static void
keyboardChar(GLFWwindow *, unsigned int codepoint)
{
    char key = (char)(codepoint & 0x7f);//unicode to ascii

    if (g_hud.KeyDown(tolower(key))) return;

    switch (key) {
        case ' ': startRender(); break;
        case 'q': g_running = 0; break;
        case 'f': fitFrame(); break;
        case '+':
        case '=': g_preTessLevel++; updateGeom(); break;
        case '-': g_preTessLevel = std::max(1, g_preTessLevel-1); updateGeom(); break;
        case 'G': dumpCamera(); break;
        case 'g': loadCamera(); setCamera(); break;
    }
}
static void
keyboard(GLFWwindow *, int key, int /* scancode */, int event, int /* mods */)
{
    if (event == GLFW_RELEASE) return;

    const float panStep = g_dolly*0.0001;
    switch (key) {
    case GLFW_KEY_UP:
        g_pan[1] -= panStep; setCamera(); break;
    case GLFW_KEY_DOWN:
        g_pan[1] += panStep; setCamera(); break;
    case GLFW_KEY_LEFT:
        g_pan[0] += panStep; setCamera(); break;
    case GLFW_KEY_RIGHT:
        g_pan[0] -= panStep; setCamera(); break;

    case GLFW_KEY_ESCAPE: g_hud.SetVisible(!g_hud.IsVisible()); display(); break;
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
callbackIntersect(int b)
{
    g_intersectKernel = b;

    startRender();
}

static void
callbackSlider(float value, int data)
{
    if (data == 0) {
        g_uvMargin = value;
        startRender();
    } else if (data == 1) {
        g_displaceScale = value;
        updateGeom();
    } else if (data == 2) {
        g_displaceFreq = value;
        startRender();
    }
}

static void
callbackCheckBox(bool checked, int button)
{
    switch (button) {
    case kHUD_CB_DISPLAY_BVH:
        g_drawBVH = checked;
        break;
    case kHUD_CB_BLOCK_FILL:
        g_blockFill = checked;
        break;
    case kHUD_CB_PRE_TESSELLATE:
        g_preTess = checked;
        updateGeom();
        break;
    case kHUD_CB_ANIMATE:
        g_animate = checked;
        updateGeom();
        break;
    case kHUD_CB_WATERTIGHT:
        g_watertight = checked;
        updateGeom();
        break;
    case kHUD_CB_CROP_UV:
        g_cropUV = checked;
        updateGeom();
        break;
    case kHUD_CB_BEZIER_CLIP:
        g_bezierClip = checked;
        updateGeom();
        break;
    case kHUD_CB_DEBUG:
        g_debug = checked;
        updateGeom();
        break;
    }
    display();
}

static void
initHUD()
{
#if GLFW_VERSION_MAJOR>=3
    // window size might not match framebuffer size on a high DPI display
    glfwGetWindowSize(g_window, &g_width, &g_height);
#endif
    g_hud.Init(g_width, g_height);

    g_hud.AddCheckBox("Show BVH (B)", g_drawBVH != 0,
                      10, 10, callbackCheckBox, kHUD_CB_DISPLAY_BVH, 'b');
    g_hud.AddCheckBox("Block Fill (K)", g_blockFill != 0,
                      10, 30, callbackCheckBox, kHUD_CB_BLOCK_FILL, 'k');

    g_hud.AddCheckBox("Watertight (C)", g_watertight != 0,
                      10, 60, callbackCheckBox, kHUD_CB_WATERTIGHT, 'c');

    g_hud.AddCheckBox("Crop UV (U)", g_cropUV != 0,
                      10, 80, callbackCheckBox, kHUD_CB_CROP_UV, 'u');

    g_hud.AddCheckBox("Bezier Clip (J)", g_bezierClip != 0,
                      10, 100, callbackCheckBox, kHUD_CB_BEZIER_CLIP, 'j');

    g_hud.AddCheckBox("Pre tessellate (T)", g_preTess != 0,
                      10, 120, callbackCheckBox, kHUD_CB_PRE_TESSELLATE, 't');
    g_hud.AddCheckBox("Animate vertices (M)", g_animate != 0,
                      10, 140, callbackCheckBox, kHUD_CB_ANIMATE, 'm');
    g_hud.AddCheckBox("Debug (D)", g_debug != 0,
                      10, 160, callbackCheckBox, kHUD_CB_DEBUG, 'd');

    g_hud.AddSlider("UV Margin", 0, 0.01, g_uvMargin,
                    10, 180, 20, false, callbackSlider, 0);
    g_hud.AddSlider("Disp scale", 0, 0.1, g_displaceScale,
                    10, 210, 20, false, callbackSlider, 1);
    g_hud.AddSlider("Disp freq", 0, 200, g_displaceFreq,
                    10, 240, 20, false, callbackSlider, 2);

    int kernel_pulldown = g_hud.AddPullDown("Intersect (I)", 400, 10, 200, callbackIntersect, 'i');
    g_hud.AddPullDownButton(kernel_pulldown, "Original", 0, g_intersectKernel == 0);
    g_hud.AddPullDownButton(kernel_pulldown, "Osd float", 1, g_intersectKernel == 1);
    g_hud.AddPullDownButton(kernel_pulldown, "Osd sse", 2, g_intersectKernel == 2);
    g_hud.AddPullDownButton(kernel_pulldown, "Osd double", 3, g_intersectKernel == 3);

    int shading_pulldown = g_hud.AddPullDown("Shading (W)", 200, 10, 250, callbackDisplayStyle, 'w');
    g_hud.AddPullDownButton(shading_pulldown, "Shaded", Scene::SHADED,
                            g_displayStyle==Scene::SHADED);
    g_hud.AddPullDownButton(shading_pulldown, "Ptex coord", Scene::PTEX_COORD,
                            g_displayStyle==Scene::PTEX_COORD);
    g_hud.AddPullDownButton(shading_pulldown, "Patch color", Scene::PATCH_TYPE,
                            g_displayStyle==Scene::PATCH_TYPE);
    g_hud.AddPullDownButton(shading_pulldown, "Clip level", Scene::CLIP_LEVEL,
                            g_displayStyle==Scene::CLIP_LEVEL);
    g_hud.AddPullDownButton(shading_pulldown, "Ambient Occlusion", Scene::AO,
                            g_displayStyle==Scene::AO);
    g_hud.AddPullDownButton(shading_pulldown, "Transparent", Scene::TRANSPARENT,
                            g_displayStyle==Scene::TRANSPARENT);

    for (int i = 1; i < 11; ++i) {
        char level[16];
        sprintf(level, "Lv. %d", i);
        g_hud.AddRadioButton(3, level, i==g_level, 10, 250+i*20, callbackLevel, i, '0'+(i%10));
    }

    int pulldown_handle = g_hud.AddPullDown("Shape (N)", -300, 10, 300, callbackModel, 'n');
    for (int i = 0; i < (int)g_defaultShapes.size(); ++i) {
        g_hud.AddPullDownButton(pulldown_handle, g_defaultShapes[i].name.c_str(), i,
                                i == g_currentShape);
    }   
}

//------------------------------------------------------------------------------

static GLuint
compileProgram(const char *vs, const char *gs, const char *fs)
{
    GLuint vertexShader = compileShader(GL_VERTEX_SHADER, vs);
    GLuint geometryShader = gs ? compileShader(GL_GEOMETRY_SHADER, gs) : 0;
    GLuint fragmentShader = compileShader(GL_FRAGMENT_SHADER, fs);
    GLuint program = glCreateProgram();
    glAttachShader(program, vertexShader);
    glAttachShader(program, fragmentShader);
    if (geometryShader) glAttachShader(program, geometryShader);
    glLinkProgram(program);
    glDeleteShader(vertexShader);
    glDeleteShader(fragmentShader);
    if (geometryShader) glDeleteShader(geometryShader);
    return program;
}

static void
initGL()
{
    glClearColor(0.1f, 0.1f, 0.1f, 1.0f);
    glEnable(GL_DEPTH_TEST);
    glDepthFunc(GL_LEQUAL);
    glCullFace(GL_BACK);
    glEnable(GL_CULL_FACE);

    glGenVertexArrays(1, &g_vao);
    glGenVertexArrays(1, &g_vaoBVH);

    glGenBuffers(1, &g_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, g_vbo);
    static float pos[] = { -1, -1,
                            1, -1,
                           -1, 1,
                            1, 1,
                            0.5, -1,
                            1, -1,
                            0.5, -0.5,
                            1, -0.5 };
    glBufferData(GL_ARRAY_BUFFER, sizeof(pos), pos, GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);

    g_programBlockFill = compileProgram(s_VS, NULL, s_FS);
    g_programSimpleFill = compileProgram(s_VS, NULL, s0_FS);
    g_programBVH = compileProgram(s_VS_BVH, s_GS_BVH, s_FS_BVH);
    g_programDebug = compileProgram(s_VS_Debug, NULL, s_FS);
}

//------------------------------------------------------------------------------
static void
idle() {

    if (g_stepIndex == 0) {
        if (g_animate) {
            ++g_frame;
            updateGeom();
        }
        return;
    }

    double fov = 45.0f;
    g_renderTime = -1.0f;

    int index = g_step*g_step - g_stepIndex;

    index = ((index>>0)&1) * (g_step*g_step>>1)
          + ((index>>1)&1) * (g_step>>1)
          + ((index>>2)&1) * (g_step*g_step>>2)
          + ((index>>3)&1) * (g_step>>2)
          + ((index>>4)&1) * (g_step*g_step>>3)
          + ((index>>5)&1) * (g_step>>3);


    g_scene.Render(g_width, g_height, fov,
                   g_image,
                   g_eye, g_lookat, g_up, g_step, index,
                   g_intersectKernel, g_uvMargin, g_cropUV, g_bezierClip,
                   g_displaceScale, g_displaceFreq);

    --g_stepIndex;

    if (g_stepIndex == 0) {
        g_renderTimer.Stop();
        g_renderTime = g_renderTimer.GetElapsed();
    }

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
        else if(!strcmp(argv[i], "-H")) {
            g_intersectKernel = 1;
        }
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

    glfwSetCharCallback(g_window, keyboardChar);
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

        if (not g_animate and g_stepIndex == 0) glfwWaitEvents();
    }

    uninitGL();
    glfwTerminate();
}

//------------------------------------------------------------------------------
