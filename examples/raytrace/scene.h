#ifndef SCENE_H
#define SCENE_H

#include "bvh_accel.h"
#include "camera.h"
#include <far/meshFactory.h>
#include <far/patchTables.h>
#include <osd/opengl.h>
#include <osd/vertex.h>

class CLTracer;

typedef OpenSubdiv::HbrMesh<OpenSubdiv::OsdVertex>     OsdHbrMesh;
typedef OpenSubdiv::HbrVertex<OpenSubdiv::OsdVertex>   OsdHbrVertex;
typedef OpenSubdiv::HbrFace<OpenSubdiv::OsdVertex>     OsdHbrFace;
typedef OpenSubdiv::HbrHalfedge<OpenSubdiv::OsdVertex> OsdHbrHalfedge;

class Scene
{
public:
    struct Config
    {
        Config() {
            intersectKernel = BVHAccel::NEW_FLOAT;
            uvMargin = 0.0f;
            cropUV = false;
            bezierClip = true;
            epsLevel = 4;
            maxLevel = 16;
            useTriangle = false;
            useRayDiffEpsilon = true;
            displaceScale = displaceFreq = 0.0f;
            conservativeTest = true;
            directBilinear = false;
            preTessLevel = 0;
            step = 1;
        }
        int intersectKernel;
        float uvMargin;
        bool cropUV;
        bool bezierClip;
        int epsLevel;
        int maxLevel;
        int step;
        bool useTriangle;
        bool useRayDiffEpsilon;
        bool conservativeTest;
        bool directBilinear;
        int preTessLevel;
        float displaceScale;
        float displaceFreq;


        std::string Dump() const;
    };

    Scene();
    ~Scene();

    void BezierConvert(float *vertices, int numVertices,
                       OpenSubdiv::FarPatchTables const *patchTables,
                       std::vector<int> const &farToHbr,
                       OsdHbrMesh *hbrMesh,
                       float displaceBound);

    void Tessellate(int level);

    void Build();
    void VBOBuild();

    void SetCamera(int width, int height, double fov,
                   std::vector<float> &image, // RGB
                   const float eye[3], const float lookat[3], const float up[3]);

    void SetConfig(Config const &config);

    void Render(int stepIndex, int step);
    void Render() { Render(0, 1); }

    void DebugTrace(float x, float y);

    void MakeReport(const char *filename);


    int GetNumPatches() const { return _mesh.numBezierPatches; }
    int GetNumTriangles() const { return _mesh.numTriangles; }

    void Shade(float rgba[4], const Intersection &isect, const Ray &ray, Context *context);

    enum ShadeMode { SHADED, PTEX_COORD, PATCH_TYPE, CLIP_LEVEL, QUADS, AO, TRANSPARENT };

    void SetShadeMode(ShadeMode mode) {
        _mode = mode;
    }

    enum BackgroundMode { GRADATION, WHITE, BLACK };

    void SetBackgroudMode(BackgroundMode mode) {
        _backgroundMode = mode;
    }
    BackgroundMode GetBackgroundMode() const {
        return _backgroundMode;
    }

    GLuint GetVBO() const { return _vbo; }
    int GetNumBVHNode() const { return (int)_accel.GetNodes().size(); }

    void SetWatertight(bool flag) { _watertight = flag; }

    size_t GetMemoryUsage() const {
        size_t mem = 0;
        if (_mesh.IsBezierMesh()) {
            mem += _mesh.bezierVertices.size() * sizeof(float);  // cp
            mem += _mesh.numBezierPatches/2; // (4bit per patch);
            //mem += _mesh.bezierBounds.size() * sizeof(float);   // bounds
        } else {
            mem += _mesh.triVertices.size() * sizeof(float); // verts
            mem += _mesh.faces.size() * sizeof(unsigned int); // indices
        }
        // bvh
        mem += _accel.GetNodes().size() * sizeof(BVHNode);
        mem += _accel.GetIndices().size() * sizeof(unsigned int);
        return mem;
    }

protected:
    void recordMetric(int id, std::ostream &out, Config const &config);

private:
    Camera _camera;
    Mesh _mesh;
    BVHAccel _accel;
    ShadeMode _mode;
    BackgroundMode _backgroundMode;
    bool _watertight;

    double _traverseTime;
    double _intersectTime;
    double _shadeTime;

    GLuint _vbo;
    int _width;
    int _height;
    float *_image;
    CLTracer *_clTracer;
};

#endif  // SCENE_H
