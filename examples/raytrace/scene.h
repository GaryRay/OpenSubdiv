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
#ifndef SCENE_H
#define SCENE_H

#include "bvh_accel.h"
#include "camera.h"
#include "texture.h"
#include <far/patchTables.h>
#include <osd/opengl.h>

class Scene
{
public:
    struct Config {
        Config() {
            intersectKernel = BVHAccel::OSD_FLOAT;
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

    void BuildBVH(int minLeafPrimitives);
    void BuildVBO();

    void SetCamera(int width, int height, double fov,
                   std::vector<float> &image, // RGB
                   const float eye[3], const float lookat[3], const float up[3]);

    void SetConfig(Config const &config);

    bool LoadEnvMap(const std::string &filename);
    void PBS(float rgba[4], const Intersection &isect, const Ray &ray,
             Context *context);
    // render, debug

    void Render(int stepIndex, int step);
    void Render() { Render(0, 1); }
    void DebugTrace(float x, float y);

    // shading style
    void EnvCol(float rgba[4], const OsdBezier::vec3f & dir);

    void Shade(float rgba[4], const Intersection &isect, const Ray &ray, Context *context);

    enum ShadeMode { SHADED, PTEX_COORD, PATCH_TYPE, HEAT_MAP, QUADS, AO, TRANSPARENT };

    void SetShadeMode(ShadeMode mode) {
        _mode = mode;
    }

    // background style

    enum BackgroundMode { GRADATION, WHITE, BLACK, ENVMAP };

    void SetBackgroudMode(BackgroundMode mode) {
        _backgroundMode = mode;
    }
    BackgroundMode GetBackgroundMode() const {
        return _backgroundMode;
    }

    // mesh, vbo

    Mesh &GetMesh() { return _mesh; }

    GLuint GetVBO() const { return _vbo; }

    int GetNumBVHNode() const { return (int)_accel.GetNodes().size(); }

    // reporting

    void RenderReport();

    size_t GetMemoryUsage() const {
        size_t mem = 0;
        // bvh
        mem += _accel.GetNodes().size() * sizeof(BVHNode);
        mem += _accel.GetIndices().size() * sizeof(unsigned int);
        return mem;
    }

private:
    Camera _camera;
    Mesh _mesh;
    BVHAccel _accel;
    ShadeMode _mode;
    BackgroundMode _backgroundMode;

    Texture _envMap;

    double _traverseTime;
    double _intersectTime;
    double _shadeTime;

    GLuint _vbo;
    int _width;
    int _height;
    float *_image;
};

#endif  // SCENE_H
