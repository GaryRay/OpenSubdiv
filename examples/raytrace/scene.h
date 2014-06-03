#ifndef SCENE_H
#define SCENE_H

#include "bvh_accel.h"
#include "mesh.h"
#include <far/patchTables.h>

class Scene
{
public:
    Scene();

    void Build(float *vertices, int numVertices,
               OpenSubdiv::FarPatchTables const *patchTables);

    void Render(int width, int height, double fov,
                std::vector<float> &image, // RGB
                std::vector<int> &count, const float eye[3],
                const float lookat[3], const float up[3],
                int step);

    int GetNumPatches() const { return _mesh.numBezierPatches; }

    void Shade(float rgba[4], const Intersection &isect, const Ray &ray);

private:
    Mesh _mesh;
    BVHAccel _accel;
};

#endif  // SCENE_H
