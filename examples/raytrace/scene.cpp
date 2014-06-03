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

unsigned int gSeed[1024][4];

inline void init_randomreal(void) {
#if _OPENMP
  assert(omp_get_max_threads() < 1024);

  for (int i = 0; i < omp_get_max_threads(); i++) {
    gSeed[i][0] = 123456789 + i;
    gSeed[i][1] = 362436069;
    gSeed[i][2] = 521288629;
    gSeed[i][3] = 88675123;
  }
#else
#endif
}

inline double randomreal(void) {
// xorshift RNG
#ifdef _OPENMP
  int tid = omp_get_thread_num();
  unsigned int x = gSeed[tid][0];
  unsigned int y = gSeed[tid][1];
  unsigned int z = gSeed[tid][2];
  unsigned int w = gSeed[tid][3];
  unsigned t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));

  gSeed[tid][0] = x;
  gSeed[tid][1] = y;
  gSeed[tid][2] = z;
  gSeed[tid][3] = w;
  return w * (1.0 / 4294967296.0);
#else
  // @fixme { don't use __thread keyword? }
  static unsigned int  x = 123456789, y = 362436069, z = 521288629,
                                 w = 88675123;
  unsigned t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  return w * (1.0 / 4294967296.0);
#endif
}


Scene::Scene()
{
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

    _mesh.numVertices = numVertices;
    _mesh.numBezierPatches = numPatches;
    _mesh.vertices = &vertices[0];
    _mesh.bezierVertices = &bezierVertices[0];
    _mesh.patchParams = &(patchParam[0]);

    assert(numPatches*16*3 == (int)bezierVertices.size());

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

  init_randomreal();

#pragma omp parallel for schedule(dynamic, 1)
  for (int y = stepIndex/step; y < height; y += step) {

    for (int x = stepIndex%step; x < width; x += step) {

        float u = 0.5;//randomreal() - 0.5;
        float v = 0.5;//randomreal() - 0.5;

      Ray ray = camera.GenerateRay(x + u + step / 2.0f, y + v + step / 2.0f);

      Intersection isect;
      bool hit = _accel.Traverse(isect, &_mesh, ray);

      if (hit) {

        float rgba[4];
        Shade(rgba, isect, ray);
        image[4 * (y * width + x) + 0] = rgba[0];
        image[4 * (y * width + x) + 1] = rgba[1];
        image[4 * (y * width + x) + 2] = rgba[2];
      } else {
          image[4 * (y * width + x) + 0] = 0.1;
          image[4 * (y * width + x) + 1] = 0.1;
          image[4 * (y * width + x) + 2] = 0.1;
      }
      image[4 * (y * width + x) + 3] = 1.0;
#if 0
      // block fill
      for (int v = 0; v < step; v++) {
          if (y+v >= height) continue;
          for (int u = 0; u < step; u++) {
              if (x+u >= width) continue;
              for (int k = 0; k < 3; k++) {
                  image[((y + v) * width * 3 + (x + u) * 3) + k] =
                      image[3 * (y * width + x) + k];
                  count[(y+v) * width + (x+u)]++;
              }
          }
      }
#endif
    }
  }

}

void
Scene::Shade(float rgba[4], const Intersection &isect, const Ray &ray)
{
    real3 I(0.4, 1, -0.4);
    I.normalize();

    Ray sray;

    real3 hitP = ray.org + isect.t * ray.dir;

    sray.org = hitP - 0.00001 * I;
    sray.dir = -1*I;

    Intersection is;
    bool hit = false; //_accel.Traverse(is, &_mesh, sray);
    if (!hit) {
        real IdotN = vdot(I, isect.normal);
        IdotN = std::max(0.2f, IdotN);
        rgba[0] = IdotN;
        rgba[1] = IdotN;
        rgba[2] = IdotN;
        rgba[3] = 1.0;
    } else {
        rgba[0] = 0.2;
        rgba[1] = 0.2;
        rgba[2] = 0.2;
        rgba[3] = 1.0;
    }
}
