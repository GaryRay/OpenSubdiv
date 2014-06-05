#ifndef SCENE_H
#define SCENE_H

#include "bvh_accel.h"
#include <far/patchTables.h>
#include <osd/opengl.h>

class Scene
{
public:
    Scene();
    ~Scene();

    void Build(float *vertices, int numVertices,
               OpenSubdiv::FarPatchTables const *patchTables);

    void Render(int width, int height, double fov,
                std::vector<float> &image, // RGB
                const float eye[3], const float lookat[3], const float up[3],
                int step, int stepIndex);

    int GetNumPatches() const { return _mesh.numBezierPatches; }

    void Shade(float rgba[4], const Intersection &isect, const Ray &ray);

    enum ShadeMode { SHADED, PTEX_COORD, PATCH_TYPE };

    void SetShadeMode(ShadeMode mode) {
        _mode = mode;
    }

    GLuint GetVBO() const { return _vbo; }
    int GetNumBVHNode() const { return (int)_accel.GetNodes().size(); }

private:
    Mesh _mesh;
    BVHAccel _accel;
    ShadeMode _mode;

    GLuint _vbo;
};

#endif  // SCENE_H
