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

#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>

Scene::Scene()
{
    static tbb::task_scheduler_init init;
}

void
Scene::Build(float *inVertices, int numVertices, OpenSubdiv::FarPatchTables const *patchTables)
{
    using namespace OpenSubdiv;

    // convert to mesh
    FarPatchTables::PatchArrayVector const &patchArrays = patchTables->GetPatchArrayVector();
    FarPatchTables::PatchParamTable const &patchParam = patchTables->GetPatchParamTable();

    // centering/normalize vertices.
    std::vector<float> vertices;
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

    int numPatches = 0;
    static std::vector<float> bezierVertices;
    bezierVertices.clear();
    // iterate patch types.
    for (FarPatchTables::PatchArrayVector::const_iterator it = patchArrays.begin();
         it != patchArrays.end(); ++it) {

        switch(it->GetDescriptor().GetType()) {
        case FarPatchTables::REGULAR:
            numPatches += convertRegular(bezierVertices, &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::BOUNDARY:
            numPatches += convertBoundary(bezierVertices, &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::CORNER:
            numPatches += convertCorner(bezierVertices, &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::GREGORY:
            numPatches += convertGregory(bezierVertices, &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::GREGORY_BOUNDARY:
            numPatches += convertBoundaryGregory(bezierVertices, &vertices[0], patchTables, *it);
            break;
        default:
            break;
        }
    }

    _mesh.numBezierPatches = numPatches;
    _mesh.bezierVertices = &bezierVertices[0];
    _mesh.patchParams = &(patchParam[0]);

    assert(numPatches*16*3 == (int)bezierVertices.size());

    BVHBuildOptions options; // Use default option

    options.minLeafPrimitives = 64;

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

class Kernel {
public:
    Kernel(int width, int stepIndex, int step, BVHAccel *accel, Mesh *mesh, Camera *camera, float *image, int mode) :
        _width(width), _stepIndex(stepIndex), _step(step), _accel(accel), _mesh(mesh), _camera(camera), _image(image), _mode(mode) {
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

            float rgba[4];
            if (hit) {
                Shade(rgba, isect, ray);
            } else {
                rgba[0] = rgba[1] = rgba[2] = 0.1f;
                rgba[3] = 1.0f;
            }
            _image[4 * (y * _width + x) + 0] = rgba[0];
            _image[4 * (y * _width + x) + 1] = rgba[1];
            _image[4 * (y * _width + x) + 2] = rgba[2];
            _image[4 * (y * _width + x) + 3] = rgba[3];
        }
        }
    }

    void Shade(float rgba[4], const Intersection &isect, const Ray &ray) const {
        real3 I = ray.dir;

        real3 color(0.8f, 0.8f, 0.8f);
        if (_mode == 1) {
            color = real3(isect.u, isect.v, 1);
        }

        real d = vdot(I, isect.normal);
        real3 reflect = I - 2 * d * isect.normal;
        real s = pow(std::max(0.0f, -vdot(ray.dir, reflect)), 32);
        rgba[0] = d * color[0] + s;
        rgba[1] = d * color[1] + s;
        rgba[2] = d * color[2] + s;
        rgba[3] = 1.0;
    }

    private:
    int _width;
      int _stepIndex;
    int _step;
      BVHAccel *_accel;
    Mesh *_mesh;
    Camera *_camera;
    float *_image;
    int _mode;
};

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
    assert(image.size() >= 3 * width * height);

#if 1
    tbb::blocked_range<int> range(0, height/step, 1);

    Kernel kernel(width, stepIndex, step, &_accel, &_mesh, &camera, &image[0], _mode);
    tbb::parallel_for(range, kernel);

#else
#pragma omp parallel for schedule(dynamic, 1)
    for (int y = stepIndex/step; y < height; y += step) {

        for (int x = stepIndex%step; x < width; x += step) {

            float u = 0.5;
            float v = 0.5;

            Ray ray = camera.GenerateRay(x + u + step / 2.0f, y + v + step / 2.0f);

            Intersection isect;
            bool hit = _accel.Traverse(isect, &_mesh, ray);

            float rgba[4];
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

void
Scene::Shade(float rgba[4], const Intersection &isect, const Ray &ray)
{
    real3 I = ray.dir;

    real3 color(0.8f, 0.8f, 0.8f);
    if (_mode == 1) {
        color = real3(isect.u, isect.v, 1);
    }

    real d = vdot(I, isect.normal);
    real3 reflect = I - 2 * d * isect.normal;
    real s = pow(std::max(0.0f, -vdot(ray.dir, reflect)), 32);
    rgba[0] = d * color[0] + s;
    rgba[1] = d * color[1] + s;
    rgba[2] = d * color[2] + s;
    rgba[3] = 1.0;
}
