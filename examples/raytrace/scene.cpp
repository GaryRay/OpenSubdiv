#include "scene.h"
#include "convert_bezier.h"
#include "camera.h"

#include <ctime>
#include <cstring>
#include <string>
#include <cfloat>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OPENSUBDIV_HAS_TBB
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#endif

static float const *getAdaptivePatchColor(OpenSubdiv::FarPatchTables::Descriptor const & desc)
{
    static float _colors[4][5][4] = {{{1.0f,  1.0f,  1.0f,  1.0f},   // regular
                                      {0.8f,  0.0f,  0.0f,  1.0f},   // boundary
                                      {0.0f,  1.0f,  0.0f,  1.0f},   // corner
                                      {1.0f,  1.0f,  0.0f,  1.0f},   // gregory
                                      {1.0f,  0.5f,  0.0f,  1.0f}},  // gregory boundary

                                     {{0.0f,  1.0f,  1.0f,  1.0f},   // regular pattern 0
                                      {0.0f,  0.5f,  1.0f,  1.0f},   // regular pattern 1
                                      {0.0f,  0.5f,  0.5f,  1.0f},   // regular pattern 2
                                      {0.5f,  0.0f,  1.0f,  1.0f},   // regular pattern 3
                                      {1.0f,  0.5f,  1.0f,  1.0f}},  // regular pattern 4
 
                                     {{0.0f,  0.0f,  0.75f, 1.0f},   // boundary pattern 0
                                      {0.0f,  0.2f,  0.75f, 1.0f},   // boundary pattern 1
                                      {0.0f,  0.4f,  0.75f, 1.0f},   // boundary pattern 2
                                      {0.0f,  0.6f,  0.75f, 1.0f},   // boundary pattern 3
                                      {0.0f,  0.8f,  0.75f, 1.0f}},  // boundary pattern 4
 
                                     {{0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 0
                                      {0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 1
                                      {0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 2
                                      {0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 3
                                      {0.25f, 0.25f, 0.25f, 1.0f}}}; // corner pattern 4

    typedef OpenSubdiv::FarPatchTables FPT;

    if (desc.GetPattern()==FPT::NON_TRANSITION) {
        return _colors[0][(int)(desc.GetType()-FPT::REGULAR)];
    } else {
        return _colors[(int)(desc.GetType()-FPT::REGULAR)+1][(int)desc.GetPattern()-1];
    }
}

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
Scene::Convert(float *inVertices, int numVertices, OpenSubdiv::FarPatchTables const *patchTables)
{
    using namespace OpenSubdiv;

    // convert to mesh
    FarPatchTables::PatchArrayVector const &patchArrays = patchTables->GetPatchArrayVector();
    FarPatchTables::PatchParamTable const &patchParam = patchTables->GetPatchParamTable();

    // centering/normalize vertices.
    std::vector<float> vertices;
    vertices.reserve(numVertices*3);
    {
        float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
        float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (int i = 0; i < numVertices; ++i) {
            float *v = inVertices + i*3;
            for (int j = 0; j < 3; ++j) {
                min[j] = std::min(min[j], v[j]);
                max[j] = std::max(max[j], v[j]);
            }
        }
        float center[3] = { (max[0]+min[0])*0.5f,
                            (max[1]+min[1])*0.5f,
                            (max[2]+min[2])*0.5f };
        float radius = std::max(std::max(max[0]-min[0], max[1]-min[1]), max[2]-min[2]);
        for (int i = 0; i < numVertices; ++i) {
            float *v = inVertices + i*3;
            vertices.push_back((v[0]-center[0])/radius);
            vertices.push_back((v[1]-center[1])/radius);
            vertices.push_back((v[2]-center[2])/radius);
        }
    }

    int numTotalPatches = 0;
    _mesh.bezierVertices.clear();
    _mesh.bezierBounds.clear();
    _mesh.colors.clear();

    // iterate patch types.
    for (FarPatchTables::PatchArrayVector::const_iterator it = patchArrays.begin();
         it != patchArrays.end(); ++it) {

        int numPatches = 0;
        FarPatchTables::Descriptor desc = it->GetDescriptor();
        switch(desc.GetType()) {
        case FarPatchTables::REGULAR:
            numPatches = convertRegular(_mesh.bezierVertices,
                                        _mesh.bezierBounds,
                                        &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::BOUNDARY:
            numPatches = convertBoundary(_mesh.bezierVertices,
                                         _mesh.bezierBounds,
                                         &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::CORNER:
            numPatches = convertCorner(_mesh.bezierVertices,
                                       _mesh.bezierBounds,
                                       &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::GREGORY:
            numPatches = convertGregory(_mesh.bezierVertices,
                                        _mesh.bezierBounds,
                                        &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::GREGORY_BOUNDARY:
            numPatches = convertBoundaryGregory(_mesh.bezierVertices,
                                                _mesh.bezierBounds,
                                                &vertices[0], patchTables, *it);
            break;
        default:
            break;
        }
        numTotalPatches += numPatches;
        // color array
        {
            const float *color = getAdaptivePatchColor(desc);
            for (int i = 0; i < numPatches; ++i) {
                _mesh.colors.push_back(color[0]);
                _mesh.colors.push_back(color[1]);
                _mesh.colors.push_back(color[2]);
            }
        }
    }

    _mesh.numBezierPatches = numTotalPatches;
    _mesh.patchParams = &(patchParam[0]);

    assert(numTotalPatches*16*3 == (int)_mesh.bezierVertices.size());
}

void
Scene::Build()
{
    BVHBuildOptions options; // Use default option

    printf("  BVH build option:\n");
    printf("    # of leaf primitives: %d\n", options.minLeafPrimitives);
    printf("    SAH binsize         : %d\n", options.binSize);

    _accel = BVHAccel();
    _accel.Build(&_mesh, options);

    BVHBuildStatistics stats = _accel.GetStatistics();

    printf("  BVH statistics:\n");
    printf("    # of leaf   nodes: %d\n", stats.numLeafNodes);
    printf("    # of branch nodes: %d\n", stats.numBranchNodes);
    printf("  Max tree depth   : %d\n", stats.maxTreeDepth);
}

void
Scene::VBOBuild()
{
    // create vbo for display
    if (_vbo == 0) glGenBuffers(1, &_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo);
    glBufferData(GL_ARRAY_BUFFER, _accel.GetNodes().size() * sizeof(BVHNode),
                 &_accel.GetNodes()[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

#ifdef OPENSUBDIV_HAS_TBB
class Kernel {
public:
    Kernel(int width, int stepIndex, int step, BVHAccel *accel, Mesh *mesh,
           Camera *camera, float *image, Scene *scene) :
        _width(width), _stepIndex(stepIndex), _step(step), _accel(accel),
        _mesh(mesh), _camera(camera), _image(image), _scene(scene) {
    }

    void operator() (tbb::blocked_range<int> const &r) const {
        for (int rr = r.begin(); rr < r.end(); ++rr) {
            int y = rr*_step + _stepIndex/_step;
            for (int x = _stepIndex%_step; x < _width; x += _step) {

                float u = 0.5;
                float v = 0.5;

                Ray ray = _camera->GenerateRay(x + u + _step / 2.0f, y + v + _step / 2.0f);

                Intersection isect;
                bool hit = _accel->Traverse(isect, _mesh, ray);

                float rgba[4] = { 0, 0, 0, 0 };
                float *d = _image + 4 * (y * _width + x);
                if (hit) {
                    _scene->Shade(rgba, isect, ray);
                } else {
                    rgba[0] = rgba[1] = rgba[2] = 0.1f;
                    rgba[3] = 1.0f;
                }
                d[0] = rgba[0];
                d[1] = rgba[1];
                d[2] = rgba[2];
                d[3] = rgba[3];
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
#endif

void
Scene::Render(int width, int height, double fov,
              std::vector<float> &image, // RGB
              const float eye[3],
              const float lookat[3], const float up[3],
              int step, int stepIndex)
{
    std::vector<int> xs;
    std::vector<int> ys;
    std::srand(unsigned(std::time(0)));
    std::random_shuffle(xs.begin(), xs.end());
    std::random_shuffle(ys.begin(), ys.end());

    Camera camera;

    double deye[3] = { eye[0], eye[1], eye[2] };
    double dlook[3] = { lookat[0], lookat[1], lookat[2] };
    double dup[3] = { up[0], up[1], up[2] };

    camera.BuildCameraFrame(deye, dlook, dup, fov, width, height);
    assert((int)image.size() >= 3 * width * height);

#ifdef OPENSUBDIV_HAS_TBB
    tbb::blocked_range<int> range(0, height/step, 1);

    Kernel kernel(width, stepIndex, step, &_accel, &_mesh, &camera, &image[0], this);
    tbb::parallel_for(range, kernel);
#else

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int y = stepIndex/step; y < height; y += step) {

        for (int x = stepIndex%step; x < width; x += step) {

            float u = 0.5;
            float v = 0.5;

            Ray ray = camera.GenerateRay(x + u + step / 2.0f, y + v + step / 2.0f);

            Intersection isect;
            bool hit = _accel.Traverse(isect, &_mesh, ray);

            float rgba[4] = { 0, 0, 0, 0};
            if (hit) {
                Shade(rgba, isect, ray);
            } else {
                rgba[0] = rgba[1] = rgba[2] = 0.1f;
                rgba[3] = 1.0f;
            }
            image[4 * (y * width + x) + 0] = rgba[0];
            image[4 * (y * width + x) + 1] = rgba[1];
            image[4 * (y * width + x) + 2] = rgba[2];
            image[4 * (y * width + x) + 3] = rgba[3];
        }
    }
#endif
}

inline float randomreal(void) {
  static unsigned int x = 123456789, y = 362436069, z = 521288629,
      w = 88675123;
  unsigned t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  return w * (1.0f / 4294967296.0f);
}

void
Scene::Shade(float rgba[4], const Intersection &isect, const Ray &ray)
{
    real3 I = ray.dir;

    real d = std::max(real(0), vdot(I, isect.normal));
    real3 color;
    if (_mode == SHADED) {
        real3 reflect = I - 2 * d * isect.normal;
        real s = pow(std::max(0.0f, -vdot(ray.dir, reflect)), 32);
        color = d * real3(0.8, 0.8, 0.8) + s * real3(1, 1, 1);
    } else if (_mode == PTEX_COORD) {
        color = d *real3(isect.u, isect.v, 1);
    } else if (_mode == PATCH_TYPE) {
        float l = isect.level * 0.05;
        color = d * (real3(&_mesh.colors[isect.patchID*3])
                     - real3(l, l, l));
        color[0] = std::max(real(0), color[0]);
        color[1] = std::max(real(0), color[1]);
        color[2] = std::max(real(0), color[2]);
    } else if (_mode == AO) {
        Intersection si;
        Ray sray;
        int numHits = 0;

        int numSamples = 16;
        for (int i = 0; i < numSamples; ++i) {
            real3 sample = real3(0.5-randomreal(), 0.5-randomreal(), 0.5-randomreal());
            sample.normalize();
            sray.dir = sample * (vdot(sample, isect.normal) > 0 ? -1 : 1);
            sray.invDir = sray.dir.neg();
            sray.org = ray.org + ray.dir * isect.t + sray.dir * 0.0001;
            numHits += _accel.Traverse(si, &_mesh, sray) ? 1 : 0;
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
            if (_accel.Traverse(si, &_mesh, sray)) {
                Shade(rgba, si, sray);
            } else {
                rgba[3] = 1.0;
                return;
            }
        } else {
            rgba[3] = 1.0;
        }
        return;
    }
    rgba[0] = color[0];
    rgba[1] = color[1];
    rgba[2] = color[2];
    rgba[3] = 1.0;
}

