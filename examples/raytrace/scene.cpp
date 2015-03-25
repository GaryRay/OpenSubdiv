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
#include <ctime>
#include <cstring>
#include <string>
#include <cfloat>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>

#include "scene.h"
#include "context.h"

#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"

#ifdef OPENSUBDIV_HAS_TBB
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#include <tbb/enumerable_thread_specific.h>
#endif

#include <osd/cpuEvalLimitContext.h>
#include <osd/cpuEvalLimitController.h>
#include "bezier/bezier.h"
#include "bezier/bezierIntersect.h"
#include "bezier/math.h"

using namespace OsdBezier;

std::string
Scene::Config::Dump() const
{
    std::stringstream ss;
    ss << "Kernel = " << intersectKernel << ", "
       << "CropUV = " << cropUV << ", "
       << "BezierClip = " << bezierClip << ", "
       << "Eps level = " << epsLevel << ", "
       << "Max level = " << maxLevel << ", "
       << "Use Triangle = " << useTriangle << ", "
       << "Use RayDiffEpsilon = " << useRayDiffEpsilon;

    return ss.str();
}

// ---------------------------------------------------------------------------
Scene::Scene() : _vbo(0)
{
#ifdef OPENSUBDIV_HAS_TBB
    static tbb::task_scheduler_init init;
#endif
}

Scene::~Scene()
{
    if (_vbo) glDeleteBuffers(1, &_vbo);
}

void
Scene::BuildBVH(int minLeafPrimitives)
{
    BVHBuildOptions options; // Use default option

    options.minLeafPrimitives = minLeafPrimitives;

    printf("  BVH build option:\n");
    printf("    # of leaf primitives: %d\n", options.minLeafPrimitives);
    printf("    SAH binsize         : %d\n", options.binSize);

    printf("  # of triangles : %ld\n", _mesh._numTriangles);
    printf("  # of bezier patches : %ld\n", _mesh._numBezierPatches);

    _accel = BVHAccel();
    _accel.Build(&_mesh, options);

    BVHBuildStatistics stats = _accel.GetStatistics();

    printf("  BVH statistics:\n");
    printf("    # of leaf   nodes: %d\n", stats.numLeafNodes);
    printf("    # of branch nodes: %d\n", stats.numBranchNodes);
    printf("  Max tree depth   : %d\n", stats.maxTreeDepth);
}

void
Scene::BuildVBO()
{
    // create vbo for display
    if (_vbo == 0) glGenBuffers(1, &_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo);
    glBufferData(GL_ARRAY_BUFFER, _accel.GetNodes().size() * sizeof(BVHNode),
                 &_accel.GetNodes()[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

#ifdef OPENSUBDIV_HAS_TBB
typedef tbb::enumerable_thread_specific<Context> EtsContext;
EtsContext etsContexts;
#else
typedef std::vector<Context> EtsContext;
EtsContext etsContexts(1);
#endif

static double GetTraverseTime() {
    double time = 0;
    for (EtsContext::const_iterator it = etsContexts.begin(); it != etsContexts.end(); ++it) {
        time += it->GetTraverseTime();
    }
    return time / etsContexts.size();
}
static double GetIntersectTime() {
    double time = 0;
    for (EtsContext::const_iterator it = etsContexts.begin(); it != etsContexts.end(); ++it) {
        time += it->GetIntersectTime();
    }
    return time / etsContexts.size();
}
static double GetShadeTime() {
    double time = 0;
    for (EtsContext::const_iterator it = etsContexts.begin(); it != etsContexts.end(); ++it) {
        time += it->GetShadeTime();
    }
    return time / etsContexts.size();
}
static void Reset() {
    for (EtsContext::iterator it = etsContexts.begin(); it != etsContexts.end(); ++it) {
        it->Reset();
    }
}

class Kernel {
public:
    Kernel(int width, int stepIndex, int step, BVHAccel *accel, Mesh *mesh,
           Camera *camera, float *image, Scene *scene) :
        _width(width), _stepIndex(stepIndex), _step(step), _accel(accel),
        _mesh(mesh), _camera(camera), _image(image), _scene(scene) {
    }

#ifdef OPENSUBDIV_HAS_TBB
    void operator() (tbb::blocked_range<int> const &r) const {
        bool useRayDiff = true;
        EtsContext::reference context = etsContexts.local();
        for (int rr = r.begin(); rr < r.end(); ++rr) {

#else
    void operator() (int begin, int end) const {
        bool useRayDiff = true;
        Context &context = etsContexts[0];
        for (int rr = begin; rr < end; ++rr) {
#endif
            int y = rr*_step + _stepIndex/_step;
            for (int x = _stepIndex%_step; x < _width; x += _step) {
                float u = 0.5;
                float v = 0.5;

                float offset = 10.0f;//0; // _step + 2.0f;
                Ray ray = _camera->GenerateRay(x + u + offset, y + v + offset);
                if(useRayDiff){
                    Ray rayO  = _camera->GenerateRay(x + offset, y + offset);
                    Ray rayDX = _camera->GenerateRay(x + 1 + offset, y + offset);
                    Ray rayDY = _camera->GenerateRay(x + offset, y + 1 + offset);
                    ray.dDdx = rayDX.dir-rayO.dir;
                    ray.dDdy = rayDY.dir-rayO.dir;
                    ray.hasDifferential = true;
                }

                Intersection isect;
                bool hit = _accel->Traverse(isect, _mesh, ray, &context);

                context.BeginShade();

                float *d = _image + 4 * (y * _width + x);
                float rgba[4] = { 0, 0, 0, 0 };
                rgba[0] = d[0];
                rgba[1] = d[1];
                rgba[2] = d[2];
                rgba[3] = d[3];
                if (hit) {
                    _scene->Shade(rgba, isect, ray, &context);
                } else {
                    if (_scene->GetBackgroundMode() == Scene::GRADATION) {
                      // Maya like gradation. Maybe helpful to check crack visually.
                      rgba[0] = 0.1f; 
                      rgba[1] = 0.1f;
                      rgba[2] = 0.4f * ((_width - x - 1)/(double)_width);
                      rgba[3] = 1.0f;

                      _scene->EnvCol(rgba, ray.dir);
                      rgba[0] = 0.5* rgba[0];
                      rgba[1] = 0.5* rgba[1];
                      rgba[2] = 0.5* rgba[2];
                      rgba[3] = 1.0f;

                    } else if (_scene->GetBackgroundMode() == Scene::WHITE) {
                      rgba[0] = 1.0f; 
                      rgba[1] = 1.0f;
                      rgba[2] = 1.0f;
                      rgba[3] = 1.0f;
                    } else {
                      rgba[0] = 0.0f; 
                      rgba[1] = 0.0f;
                      rgba[2] = 0.0f;
                      rgba[3] = 1.0f;
                    }
                }
                d[0] = rgba[0];
                d[1] = rgba[1];
                d[2] = rgba[2];
                d[3] = rgba[3];

                context.EndShade();
            }
        }
    }

private:
    int _width;
    int _stepIndex;
    int _step;
    BVHAccel *_accel;
    Mesh *_mesh;
    Camera *_camera;
    float *_image;
    Scene *_scene;
};

void
Scene::SetCamera(int width, int height, double fov,
                 std::vector<float> &image, // RGBA
                 const float eye[3], const float lookat[3], const float up[3])
{
    _width = width;
    _height = height;
    _image = &image[0];

    double deye[3] = { eye[0], eye[1], eye[2] };
    double dlook[3] = { lookat[0], lookat[1], lookat[2] };
    double dup[3] = { up[0], up[1], up[2] };

    _camera.BuildCameraFrame(deye, dlook, dup, fov, width, height);
    assert((int)image.size() >= 3 * width * height);
}

void
Scene::SetConfig(Config const &config)
{
    _accel.SetIntersectKernel(config.intersectKernel);
    _accel.SetUVMargin(config.uvMargin);
    _accel.SetCropUV(config.cropUV);
    _accel.SetBezierClip(config.bezierClip);
    _accel.SetDisplacement(config.displaceScale, config.displaceFreq);

    static const double EPS_FROM_LEVEL[] = {
        1e-3,1e-3,1e-3,1e-3,1e-4,1e-5,1e-6,1e-7,1e-8,1e-9,
        1e-10,1e-11,1e-12,1e-13,1e-14,1e-15,1e-16,1e-17,1e-18,1e-19
    };

    _accel.SetEpsilon (EPS_FROM_LEVEL[config.epsLevel]);
    _accel.SetMaxLevel(config.maxLevel*2);
    _accel.SetUseTriangle(config.useTriangle);
    _accel.SetUseRayDiffEpsilon(config.useRayDiffEpsilon);
    _accel.SetConservativeTest(config.conservativeTest);
    _accel.SetDirectBilinear(config.directBilinear);
}

void
Scene::Render(int stepIndex, int step)
{
#ifdef OPENSUBDIV_HAS_TBB
    tbb::blocked_range<int> range(0, _height/step, 1);

    Kernel kernel(_width, stepIndex, step, &_accel, &_mesh, &_camera, _image, this);
    Reset();
    tbb::parallel_for(range, kernel);
#else
    Kernel kernel(_width, stepIndex, step, &_accel, &_mesh, &_camera, _image, this);
    kernel(0, _height/step);
#endif

    _traverseTime += GetTraverseTime() * 1000;
    _intersectTime += GetIntersectTime() * 1000;
    _shadeTime += GetShadeTime() * 1000;
}

void
Scene::DebugTrace(float x, float y)
{
    printf("------------------------------------------\n");
    printf("Debug Trace at pixel %f, %f\n", x, y);

    float u = 0.5;
    float v = 0.5;

    Ray ray = _camera.GenerateRay(x + u, y + v);

    Intersection isect;
    bool hit = _accel.Traverse(isect, &_mesh, ray, NULL);

    float rgba[4] = { 0, 0, 0, 0};
    if (hit) {
        Shade(rgba, isect, ray, NULL);
    }
    printf("%f %f %f %f\n", rgba[0], rgba[1], rgba[2], rgba[3]);
}

void
Scene::RenderReport()
{
    int iteration = 10;

    _traverseTime = 0;
    _intersectTime = 0;
    _shadeTime = 0;

    Stopwatch s;
    s.Start();
    for (int i = 0; i < iteration; ++i) {
        Render();
    }
    s.Stop();
    float renderTime = s.GetElapsed();

    std::cout << "["
              << "Traverse = " << _traverseTime/1000.0/iteration << ", "
              << "Intersect = " << _intersectTime/1000.0/iteration << ", "
              << "Shade = " << _shadeTime/1000.0/iteration << ", "
              << "Total = " << renderTime/iteration << "], \n";
}

inline float
randomreal(void) {
    static unsigned int x = 123456789, y = 362436069, z = 521288629,
        w = 88675123;
    unsigned t = x ^ (x << 11);
    x = y;
    y = z;
    z = w;
    w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    return w * (1.0f / 4294967296.0f);
}

// Simple heat map coloring
inline void
ShadeHeatmap(float col[3], float val, float maxVal) {
    float blue[3]; blue[0] = 0.0; blue[1] = 0.0; blue[2] = 1.0;
    //vector3 green(0.0, 1.0, 0.0);
    float red[3]; red[0] = 1.0; red[1] = 0.0; red[2] = 0.0;
    //vector3 red(1.0, 0.0, 0.0);

    // 0 -> blue, 50 -> (blue + red)/2, 100 -> red
    if (val < 0.0) val = 0.0;
    if (val > maxVal) val = maxVal;
    float t = val / maxVal; // => [0, 1]

    col[0] = (1.0 - t) * blue[0] + t * red[0];
    col[1] = (1.0 - t) * blue[1] + t * red[1];
    col[2] = (1.0 - t) * blue[2] + t * red[2];
}

void
Scene::EnvCol(float rgba[4], const OsdBezier::vec3f & dir)
{
    float d[3];
    d[0] = dir[0];
    d[1] = dir[1];
    d[2] = dir[2];

    if (_envMap.IsValid()) {
        LongLatMapSampler::Sample(rgba, d, &_envMap);
        rgba[0] *= 3.14;
        rgba[1] *= 3.14;
        rgba[2] *= 3.14;
    } else {
        rgba[0] = 0.0;
        rgba[1] = 0.0;
        rgba[2] = 0.0;
        rgba[3] = 1.0;
    }
}

// apply various shadings
//
void
Scene::Shade(float rgba[4], const Intersection &isect, const Ray &ray, Context *context)
{
    OsdBezier::vec3f I = ray.dir;

    // No zero in d to distinguish crack pixel color(dark background color)
    float d = std::max(float(0.2), dot(I, isect.normal));
    float s = 0;//pow(std::max(0.0f, -dot(I - 2*d*isect.normal, I)), 64);

    OsdBezier::vec3f color;
    if (_mode == SHADED) {
        //color = d * OsdBezier::vec3f(0.8, 0.8, 0.8) + s * OsdBezier::vec3f(1, 1, 1);
        //color = ray.org + ray.dir * isect.t;
        //color[2] = color[2]  * 10;
        //color = isect.normal * 0.5 + real3(0.5, 0.5, 0.5);
        float c[4] = {0, 0, 0, 0};
        PBS(c, isect, ray, context);
        color = vec3f(c[0], c[1], c[2]);
    } else if (_mode == PTEX_COORD) {
        color = d * OsdBezier::vec3f(isect.u, isect.v, 1) + s * OsdBezier::vec3f(1.0f);
        // float m = 0.01;
        // float d = 1-std::max(std::max(fabs(m-isect.u), fabs((1-m)-isect.u)),
        //                      std::max(fabs(m-isect.v), fabs((1-m)-isect.v)));
        // d = d * dot(ray.dir, isect.normal);
        // if (d > 0.1f) d = 1.0f;
        // color = OsdBezier::vec3f(d);
    } else if (_mode == PATCH_TYPE) {
        float l = isect.level * 0.05;
        color = d * (OsdBezier::vec3f(&_mesh._colors[isect.patchID*3])
                     - OsdBezier::vec3f(l, l, l));
        color[0] = std::max(0.0f, color[0]);
        color[1] = std::max(0.0f, color[1]);
        color[2] = std::max(0.0f, color[2]);
    } else if (_mode == HEAT_MAP) {
        float col[3];
        ShadeHeatmap(col, isect.clipLevel, isect.maxLevel);
        //printf("col = %f, %f, %f(lv:%d, mlv: %d)\n", col[0], col[1], col[2], isect.level, isect.maxLevel);
        color[0] = col[0];
        color[1] = col[1];
        color[2] = col[2];
    } else if (_mode == QUADS) {
        color[0] = d*((((isect.quadHash>>0)&0xff)/255.0)*0.5 + 0.5);
        color[1] = d*((((isect.quadHash>>8)&0xff)/255.0)*0.5 + 0.5);
        color[2] = d*((((isect.quadHash>>16)&0xff)/255.0)*0.5 + 0.5);
    } else if (_mode == AO) {
        Intersection si;
        Ray sray;
        int numHits = 0;

        int numSamples = 16;
        for (int i = 0; i < numSamples; ++i) {
            OsdBezier::vec3f sample(0.5-randomreal(), 0.5-randomreal(), 0.5-randomreal());
            sample.normalize();
            sray.dir = sample * (dot(sample, isect.normal) > 0 ? -1 : 1);
            sray.invDir = sray.dir.neg();
            sray.org = ray.org + ray.dir * isect.t + sray.dir * 0.0001;
            numHits += _accel.Traverse(si, &_mesh, sray, context) ? 1 : 0;
        }

        color[0] = color[1] = color[2] = d * (1.0-numHits/float(numSamples));
    } else if (_mode == TRANSPARENT) {
        float alpha = 0.25 * (1.0 - rgba[3]);
        rgba[0] += d * alpha;
        rgba[1] += d * alpha;
        rgba[2] += d * alpha;
        rgba[3] += alpha;

        if (alpha < 0.9) {
            Intersection si;
            Ray sray;
            sray.dir = ray.dir;
            sray.invDir = ray.invDir;
            sray.org = ray.org + ray.dir * (isect.t + 0.0001);
            if (_accel.Traverse(si, &_mesh, sray, context)) {
                Shade(rgba, si, sray, context);
            } else {
                rgba[3] = 1.0;
                return;
            }
        } else {
            rgba[3] = 1.0;
        }
        return;
    }
    rgba[0] += color[0];
    rgba[1] += color[1];
    rgba[2] += color[2];
    rgba[3] += 1.0;
}

bool
Scene::LoadEnvMap(const std::string &filename)
{
  int x, y, n;
  float *data = stbi_loadf(filename.c_str(), &x, &y, &n, 0);
  float gamma = 1.0f;

  Texture::Coordinate coord_type = Texture::COORDINATE_LONGLAT;

  if (data) {
    assert(n == 3); // RGB
    printf("envmap [%s] %d x %d\n", filename.c_str(), x, y);

    _envMap.Set(reinterpret_cast<const unsigned char *>(data), x, y, 3,
                Texture::FORMAT_FLOAT, gamma, coord_type);
  } else {
    printf("Failed to load envmap [%s]\n", filename.c_str());
  }

  return true;
}

float vavg(vec3f x)
{
  return (x[0] + x[1] + x[2]) / 3;
}


vec3f vclamp01(vec3f x)
{
  vec3f ret;
  ret[0] = x[0];
  ret[1] = x[1];
  ret[2] = x[2];

  return ret;
}

vec3f reflect(const vec3f &in, const vec3f &n) {
  float d = dot(in, n);
  return in - n * (2.0 * d);
}

vec3f refract(bool &tir, const vec3f &in, const vec3f &n, float eta) {
  vec3f ret;
  vec3f N;
  double e = eta;
  double cos1 = dot(in, n);
  if (cos1 < 0.0) { // entering
    N = n;
  } else { // outgoing
    cos1 = -cos1;
    e = 1.0f / eta;
    N = n.neg();
  }

  double k = 1.0 - (e * e) * (1.0 - cos1 * cos1);
  if (k <= 0.0) {
    // Toral internal reflection.
    ret = reflect(in, n);
    tir = true;
    ret.normalize();
    return ret;
  }

  k = -e * cos1 - sqrt(k);

  tir = false;
  ret = k * N + e * in;
  ret.normalize();

  return ret;
}

void fresnel_factor(vec3f &refl, vec3f& refr, float &kr, float &kt,
             const vec3f &in, const vec3f &n, float eta) {
  float d = dot(in, n);

  refl = reflect(in, n);

  bool tir;
  refr = refract(tir, in, n, eta);

  if (tir) {
    kr = 1.0;
    kt = 0.0;
    return;
  }

  float cos_r = dot(refl, n);
  float cos_t = -dot(refr, n);

  float rp = (cos_t - eta * cos_r) / (cos_t + eta * cos_r);
  float rs = (cos_r - eta * cos_t) / (cos_r + eta * cos_t);
  kr = (rp * rp + rs * rs) * 0.5f;

  if (kr < 0.0f)
    kr = 0.0f;
  if (kr > 1.0f)
    kr = 1.0f;

  kt = 1.0f - kr;
}

void GenerateBasis(vec3f &tangent, vec3f &binormal,
                          const vec3f &normal) {
  // Find the minor axis of the vector
  int i;
  int index = -1;
  double minval = 1.0e+6;
  double val = 0;

  for (int i = 0; i < 3; i++) {
    val = fabsf(normal[i]);
    if (val < minval) {
      minval = val;
      index = i;
    }
  }

  if (index == 0) {

    tangent[0] = 0.0;
    tangent[1] = -normal[2];
    tangent[2] = normal[1];
    tangent.normalize();

    binormal = cross(tangent, normal);
    binormal.normalize();

  } else if (index == 1) {

    tangent[0] = -normal[2];
    tangent[1] = 0.0;
    tangent[2] = normal[0];
    tangent.normalize();

    binormal = cross(tangent, normal);
    binormal.normalize();

  } else {

    tangent[0] = -normal[1];
    tangent[1] = normal[0];
    tangent[2] = 0.0;
    tangent.normalize();

    binormal = cross(tangent, normal);
    binormal.normalize();
  }
}

double SampleDiffuseIS(vec3f &dir, const vec3f &normal) {
  vec3f tangent, binormal;

  GenerateBasis(tangent, binormal, normal);

  double theta = acos(sqrt(1.0 - randomreal()));
  double phi = 2.0 * M_PI * randomreal();

  double cosTheta = cos(theta);

  /* D = T*cos(phi)*sin(theta) + B*sin(phi)*sin(theta) + N*cos(theta) */
  double cos_theta = cos(theta);
  vec3f T = tangent * cos(phi) * sin(theta);
  vec3f B = binormal * sin(phi) * sin(theta);
  vec3f N = normal * (cos_theta);

  dir = T + B + N;

  return cos_theta; // PDF = weight
}

// (Modified) Ward glossy BRDF
// http://www.graphics.cornell.edu/~bjw/wardnotes.pdf
// Some from OpenShadingLanguage
void WardBRDF(
    vec3f* omega_in, // output
    float* pdf,     // output
    float* weight,  // output
    float ax,
    float ay,
    vec3f  omega_out,
    vec3f  normal,               // shading normal
    vec3f  geometric_normal)     // geometric normal
{
    float cosNO = dot(normal, omega_out);

    (*pdf) = 0.0f;
    (*weight) = 0.0f;
    (*omega_in)[0] = 0.0f;
    (*omega_in)[1] = 0.0f;
    (*omega_in)[2] = 0.0f;

    if (cosNO > 0.0f) {
        // @todo { Supply tangent vector for true aniso-brdf. }
        vec3f tangent, binormal;

        GenerateBasis(tangent, binormal, normal);

        float randu = randomreal();
        float randv = randomreal();

        float alphaRatio = ay / ax;
        float cosPhi, sinPhi;

        if (randu < 0.25f) {
            float val = 4 * randu;
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            cosPhi = 1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = tanPhi * cosPhi;
        } else if (randu < 0.5) {
            float val = 1 - 4 * (0.5f - randu);
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            // phi = (float) M_PI - phi;
            cosPhi = -1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = -tanPhi * cosPhi;
        } else if (randu < 0.75f) {
            float val = 4 * (randu - 0.5f);
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            //phi = (float) M_PI + phi;
            cosPhi = -1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = tanPhi * cosPhi;
        } else {
            float val = 1 - 4 * (1 - randu);
            float tanPhi = alphaRatio * tanf((float) M_PI_2 * val);
            // phi = 2 * (float) M_PI - phi;
            cosPhi = 1 / sqrtf(1 + tanPhi * tanPhi);
            sinPhi = -tanPhi * cosPhi;
        }

        // eq. 6
        // we take advantage of cos(atan(x)) == 1/sqrt(1+x^2)
        //                  and sin(atan(x)) == x/sqrt(1+x^2)
        float thetaDenom = (cosPhi * cosPhi) / (ax * ax) + (sinPhi * sinPhi) / (ay * ay);
        float tanTheta2 = -logf(1 - randv) / thetaDenom;
        float cosTheta  = 1 / sqrtf(1 + tanTheta2);
        float sinTheta  = cosTheta * sqrtf(tanTheta2);

        vec3f h; // already normalized becaused expressed from spherical coordinates
        h[0] = sinTheta * cosPhi;
        h[1] = sinTheta * sinPhi;
        h[1] = cosTheta;
        // compute terms that are easier in local space
        float dotx = h[0] / ax;
        float doty = h[1] / ay;
        float dotn = h[2];
        // transform to world space
        h = h[0] * tangent + h[1] * binormal + h[2] * normal;
        // generate the final sample
        float oh = dot(h, omega_out);
        //omega_in->x = 2 * oh * h.x - omega_out.x;
        //omega_in->y = 2 * oh * h.y - omega_out.y;
        //omega_in->z = 2 * oh * h.z - omega_out.z;
        (*omega_in)[0] = omega_out[0] - 2 * oh * h[0];
        (*omega_in)[1] = omega_out[1] - 2 * oh * h[1];
        (*omega_in)[1] = omega_out[2] - 2 * oh * h[2];

        float ng_dot_wi = dot(geometric_normal, (*omega_in));
        if (ng_dot_wi > 0.0f) {
            float cosNI = dot(normal, (*omega_in));
            if (cosNI > 0.0f) {
                // eq. 9
                float exp_arg = (dotx * dotx + doty * doty) / (dotn * dotn);
                float denom = 4 * (float) M_PI * ax * ay * oh * dotn * dotn * dotn;
                (*pdf) = expf(-exp_arg) / denom;
                // compiler will reuse expressions already computed
                denom = (4 * (float) M_PI * ax * ay * sqrtf(cosNO * cosNI));
                float power = cosNI * expf(-exp_arg) / denom;
                (*weight) = power;
                //domega_in_dx = (2 * m_N.dot(domega_out_dx)) * m_N - domega_out_dx;
                //domega_in_dy = (2 * m_N.dot(domega_out_dy)) * m_N - domega_out_dy;
                // Since there is some blur to this reflection, make the
                // derivatives a bit bigger. In theory this varies with the
                // roughness but the exact relationship is complex and
                // requires more ops than are practical.
                //domega_in_dx *= 10;
                //domega_in_dy *= 10;
            }
        }
    }
}

// Physically-based shader
void
Scene::PBS(float rgba[4], const Intersection &isect,
           const Ray &ray, Context *context)
{
    //printf("depth = %d\n", ray.depth);
    if (ray.depth > 3) {
        rgba[0] = 0.01;
        rgba[1] = 0.01;
        rgba[2] = 0.01;
        rgba[3] = 1.0;
        return;
    }

    // Currently available: diffuse + reflection(+glossy reflection)

    //int matID = isect.matID;
    int matID = 0;
    //const Material& mat = scene.GetMaterial(matID);

  // Preserve Energy conservation for each channel.
#if 0
    vec3f diffuse = mat.diffuse;
    vec3f reflection = mat.reflection;
    vec3f refraction = mat.refraction;
    float reflectionGlossiness = mat.reflection_glossiness;
    float refractionGlossiness = mat.refraction_glossiness;
    bool  fresnel = mat.fresnel;
    float ior = mat.ior;
#endif
    vec3f diffuse = vec3f(0.5f);
    vec3f reflection = vec3f(0.1f);
    vec3f refraction = vec3f(0.0f);
    float reflectionGlossiness = 1.0f;
    float refractionGlossiness = 1.0f;
    bool  fresnel = false;
    float ior = 0.0f;

    vec3f in, n;
    in[0] = ray.dir[0];
    in[1] = ray.dir[1];
    in[2] = ray.dir[2];
    in.normalize();

    n[0] = -isect.normal[0];
    n[1] = -isect.normal[1];
    n[2] = -isect.normal[2];
    n.normalize();

    float eta = 1.0 / ior;
    vec3f ns;

    // ks wins, kt next, kd weaks.
    vec3f one(1.0, 1.0, 1.0);
    vec3f ksRGB0 = reflection;
    vec3f ktRGB0 = refraction;
    vec3f ksRGB = ksRGB0;
    vec3f ktRGB = vclamp01((one - ksRGB0) * ktRGB0);
    vec3f kdRGB = vclamp01((one - ksRGB - ktRGB) * diffuse);

    if (fresnel) { // adjust ks and kd energy with fresnel factor.
        vec3f ns = n;
        float IdotN = dot(in, ns);
        if (IdotN < 0.0) { // outside -> inside
        } else {
            eta = ior;
            ns = n.neg();
        }
        
        vec3f sDir, tDir;
        float fresnelKr, fresnelKt;
        fresnel_factor(sDir, tDir, fresnelKr, fresnelKt, in, ns, eta);
    // sDir and tDir not used.

        ksRGB = fresnelKr * ksRGB0;
        ktRGB = fresnelKt * ktRGB0;
        kdRGB = vclamp01((one - ksRGB - ktRGB) * diffuse);
    }

    //  printf("diff = %f, %f, %f\n", diffuse[0], diffuse[1], diffuse[2]);
    //  printf("ks = %f, %f, %f\n", ksRGB[0], ksRGB[1], ksRGB[2]);
    //printf("kd = %f, %f, %f\n", kdRGB[0], kdRGB[1], kdRGB[2]);

    float ks = vavg(ksRGB); ks = std::min(1.0f, std::max(0.0f, ks));
    float kt = vavg(ktRGB); kt = std::min(1.0f, std::max(0.0f, kt));
    float kd = vavg(kdRGB); kd = std::min(1.0f, std::max(0.0f, kd));

  vec3f kdRet(0.0, 0.0, 0.0);
  vec3f ktRet(0.0, 0.0, 0.0);
  vec3f ksRet(0.0, 0.0, 0.0);
  if (kd > 0.0) {

    vec3f newDir;
    double pdf = SampleDiffuseIS(newDir, n);

    Ray diffuseRay;
    diffuseRay.dir = newDir;
    diffuseRay.invDir = newDir.neg();
    diffuseRay.org = isect.position + 0.00001f * newDir;
    //diffuseRay.org = ray.org + ray.dir * isect.t + sray.dir * 0.0001;
    diffuseRay.depth = ray.depth + 1;

    Intersection diffuseIsect;
    diffuseIsect.t = 1.0e+30;
    bool hit = _accel.Traverse(diffuseIsect, &_mesh, diffuseRay, context);
    //    bool hit = TraceRay(diffuseIsect, scene, diffuseRay);
    if (hit) {
      float diffuseRGBA[4];
      PBS(diffuseRGBA, diffuseIsect, diffuseRay, context);

      kdRet[0] = kdRGB[0] * diffuseRGBA[0];
      kdRet[1] = kdRGB[1] * diffuseRGBA[1];
      kdRet[2] = kdRGB[2] * diffuseRGBA[2];
      kdRet[3] = 1.0; // fixme
    } else {
      // env light
      float rgba[4];
      EnvCol(rgba, newDir);

      float ndotl = dot(n, newDir);
      ndotl = std::max(0.1f, ndotl);

      kdRet[0] = ndotl * kdRGB[0] * rgba[0];
      kdRet[1] = ndotl * kdRGB[1] * rgba[1];
      kdRet[2] = ndotl * kdRGB[2] * rgba[2];
    }
  }

  // reflection
  if (ks > 0.0) {
    vec3f r;

    float weight = 1.0;

    if (reflectionGlossiness < 1.0) {
      // glossy reflection. WardBRDF
      
      float pdf;

      // larget = sharper.
      float ax = 1.0f * (1.0f - reflectionGlossiness); // isotropic
      float ay = 1.0f * (1.0f - reflectionGlossiness);

      vec3f wi = in.neg();

      WardBRDF(&r, &pdf, &weight, ax, ay, wi, n, n);
      //printf("w = %f, r = %f, %f, %f, n = %f, %f, %f\n",
      //  weight, r[0], r[1], r[2], n[0], n[1], n[2]);

      // HACK
      r = r.neg();
      weight = 1.0;

    } else {
      // perfect specular.
      r = reflect(in, n);
    }

    float dir[3];
    dir[0] = r[0];
    dir[1] = r[1];
    dir[2] = r[2];

    float rmag = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

    if ((weight > 0.0) && (rmag > 0.1)) {

      Ray reflRay;
      reflRay.dir = r;
      reflRay.org = isect.position + 0.01f * r;
      reflRay.depth = ray.depth + 1;

      Intersection reflIsect;
      reflIsect.t = 1.0e+30;
      bool hit = _accel.Traverse(reflIsect, &_mesh, reflRay, context);
      //bool hit = TraceRay(reflIsect, scene, reflRay, context);

      if (hit) {

        float reflRGBA[4];
        PBS(reflRGBA, reflIsect, reflRay, context);

        ksRet[0] = ksRGB[0] * reflRGBA[0];
        ksRet[1] = ksRGB[1] * reflRGBA[1];
        ksRet[2] = ksRGB[2] * reflRGBA[2];
        ksRet[3] = 1.0; // fixme

      } else {
        // env light
        float rgba[4];
        EnvCol(rgba, r);

        ksRet[0] = ksRGB[0] * rgba[0];
        ksRet[1] = ksRGB[1] * rgba[1];
        ksRet[2] = ksRGB[2] * rgba[2];
      }

    } else {
      // ???
      ksRet[0] = 0.0;
      ksRet[1] = 0.0;
      ksRet[2] = 0.0;
    }
  }

  // refraction
  if (kt > 0.0) {

    // Simple Ward refraction.
    // @todo { GGX Glossy transmission. }
    vec3f r;

    float weight = 1.0;

    if (refractionGlossiness < 1.0) {
      // glossy reflection. WardBRDF
      
      float pdf;

      // larget = sharper.
      float ax = 1.0f * (1.0f - refractionGlossiness); // isotropic
      float ay = 1.0f * (1.0f - refractionGlossiness);

      vec3f wi = in.neg();

      WardBRDF(&r, &pdf, &weight, ax, ay, wi, n, n);
      //printf("w = %f, r = %f, %f, %f, n = %f, %f, %f\n",
      //  weight, r[0], r[1], r[2], n[0], n[1], n[2]);

      // HACK
      r = r.neg();
      weight = 1.0;

    } else {
      // perfect transmission.
      bool tir = false;
      r = refract(tir, in, n, eta);
    }

    float dir[3];
    dir[0] = r[0];
    dir[1] = r[1];
    dir[2] = r[2];

    float rmag = r[0] * r[0] + r[1] * r[1] + r[2] * r[2];

    if ((weight > 0.0) && (rmag > 0.1)) {

      Ray refrRay;
      refrRay.dir = r;
      refrRay.org = isect.position + 0.01f * r;
      refrRay.depth = ray.depth + 1;

      Intersection refrIsect;
      refrIsect.t = 1.0e+30;
      bool hit = _accel.Traverse(refrIsect, &_mesh, refrRay, context);
      //bool hit = TraceRay(refrIsect, scene, refrRay);

      if (hit) {

        float refrRGBA[4];
        PBS(refrRGBA, refrIsect, refrRay, context);

        ktRet[0] = ktRGB[0] * refrRGBA[0];
        ktRet[1] = ktRGB[1] * refrRGBA[1];
        ktRet[2] = ktRGB[2] * refrRGBA[2];
        ktRet[3] = 1.0; // fixme

      } else {
        // env light
        float rgba[4];
        EnvCol(rgba, r);

        ktRet[0] = ktRGB[0] * rgba[0];
        ktRet[1] = ktRGB[1] * rgba[1];
        ktRet[2] = ktRGB[2] * rgba[2];
      }

    } else {
      // ???
      ktRet[0] = 0.0;
      ktRet[1] = 0.0;
      ktRet[2] = 0.0;
    }
  }

  // @fixme.
  rgba[0] = 0.5*3.14 * (kdRet[0] + ksRet[0] + ktRet[0]);
  rgba[1] = 0.5*3.14 * (kdRet[1] + ksRet[1] + ktRet[1]);
  rgba[2] = 0.5*3.14 * (kdRet[2] + ksRet[2] + ktRet[2]);
  rgba[3] = 1.0;
}
