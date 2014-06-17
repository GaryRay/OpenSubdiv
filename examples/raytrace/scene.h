#ifndef SCENE_H
#define SCENE_H

#include "bvh_accel.h"
#include <far/meshFactory.h>
#include <far/patchTables.h>
#include <osd/opengl.h>
#include <osd/vertex.h>

typedef OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>     OsdHbrMesh;
typedef OpenSubdiv::HbrVertex<OpenSubdiv::OsdVertex>   OsdHbrVertex;
typedef OpenSubdiv::HbrFace<OpenSubdiv::OsdVertex>     OsdHbrFace;
typedef OpenSubdiv::HbrHalfedge<OpenSubdiv::OsdVertex> OsdHbrHalfedge;

class Scene
{
public:
    Scene();
    ~Scene();

    void BezierConvert(float *vertices, int numVertices,
                       OpenSubdiv::FarPatchTables const *patchTables,
                       std::vector<int> const &farToHbr,
                       //std::vector<int> const &vertexParentIDs,
                       OsdHbrMesh *hbrMesh);

    void Tessellate(int level);

    void Build();
    void VBOBuild();

    void Render(int width, int height, double fov,
                std::vector<float> &image, // RGB
                const float eye[3], const float lookat[3], const float up[3],
                int step, int stepIndex, int intersectKernel, float uvMargin, bool cropUV);

    int GetNumPatches() const { return _mesh.numBezierPatches; }
    int GetNumTriangles() const { return _mesh.numTriangles; }

    void Shade(float rgba[4], const Intersection &isect, const Ray &ray);

    enum ShadeMode { SHADED, PTEX_COORD, PATCH_TYPE, CLIP_LEVEL, AO, TRANSPARENT };

    void SetShadeMode(ShadeMode mode) {
        _mode = mode;
    }

    GLuint GetVBO() const { return _vbo; }
    int GetNumBVHNode() const { return (int)_accel.GetNodes().size(); }

    void SetWatertight(bool flag) { _watertight = flag; }

    size_t GetMemoryUsage() const {
        size_t mem = 0;
        if (_mesh.IsBezierMesh()) {
            mem += _mesh.bezierVertices.size() * sizeof(float);  // cp
            mem += _mesh.bezierBounds.size() * sizeof(float);   // bounds
        } else {
            mem += _mesh.triVertices.size() * 3 * sizeof(float); // verts
            mem += _mesh.faces.size() * sizeof(unsigned int); // indices
        }
        return mem;
    }


private:
    Mesh _mesh;
    BVHAccel _accel;
    ShadeMode _mode;
    bool _watertight;

    GLuint _vbo;
};

#endif  // SCENE_H
