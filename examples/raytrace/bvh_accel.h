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
#ifndef __BVH_ACCEL_H__
#define __BVH_ACCEL_H__

#include <vector>
#include "ray.h"
#include "mesh.h"
#include "context.h"
#include "bezier/math.h"

typedef float real;

class BVHNode {
public:
    typedef OsdBezier::vec3f real3;

    BVHNode() {};
    ~BVHNode() {};

    real3 bmin, bmax;
    int flag; // 1 = leaf node, 0 = branch node
    int axis;
    // leaf
    //   data[0] = npoints
    //   data[1] = index
    //
    // branch
    //   data[0] = child[0]
    //   data[1] = child[1]
    unsigned int data[2];
};

///< BVH build option.
struct BVHBuildOptions {
    bool debugPrint;
    real costTaabb;
    int minLeafPrimitives;
    int maxTreeDepth;
    int binSize;

    // Set default value: Taabb = 0.2
    BVHBuildOptions()
        : costTaabb(0.2), minLeafPrimitives(16), maxTreeDepth(256), binSize(64) {}
};

///< BVH build statistics.
struct BVHBuildStatistics {
    int maxTreeDepth;
    int numLeafNodes;
    int numBranchNodes;

    // Set default value: Taabb = 0.2
    BVHBuildStatistics() : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

typedef struct {
    typedef OsdBezier::vec3f real3;

    float t;
    float u;
    float v;
    unsigned int faceID;

    // patch info
    unsigned int patchID;
    unsigned int level;
    unsigned int clipLevel;
    unsigned int quadHash;

    // for SGA tech brief
    float        eps;
    unsigned int maxLevel;

    unsigned int f0;
    unsigned int f1;
    unsigned int f2;

    real3 position;
    real3 geometricNormal;
    real3 normal;
    real3 tangent;
    real3 binormal;
    float texcoord[2];
} Intersection;

class BVHAccel {
public:
    BVHAccel(float uvMargin=0.1f) : _intersectKernel(OSD_FLOAT),
                                    _uvMargin(uvMargin), _cropUV(true),
                                    _bezierClip(true), _displaceScale(0), _displaceFreq(100) ,
                                    _epsilon(1e-4), _maxLevel(10), _useTriangle(false),
                                    _useRayDiffEpsilon(false), _conservativeTest(false),
                                    _directBilinear(false) {}
    ~BVHAccel() {};

    ///< Build BVH for input mesh.
    bool Build(const Mesh *mesh, const BVHBuildOptions &options);

    ///< Get statistics of built BVH tree. Valid after Build()
    BVHBuildStatistics GetStatistics() const { return _stats; }

    ///< Traverse into BVH along ray and find closest hit point if found
    bool Traverse(Intersection &isect, const Mesh *mesh, Ray &ray, Context *context);

    enum { OSD_FLOAT, OSD_SSE, OSD_DOUBLE } IntersectKernel;
    void SetIntersectKernel(int k) {_intersectKernel = k; }
    void SetUVMargin(float margin) { _uvMargin = margin; }
    void SetCropUV(bool flag) {_cropUV = flag;}
    void SetBezierClip(bool flag) {_bezierClip = flag;}
    void SetDisplacement(float scale, float freq) { _displaceScale = scale; _displaceFreq = freq; }

    void SetEpsilon(double eps){_epsilon=eps;}
    void SetMaxLevel(int level){_maxLevel=level;}
    void SetUseTriangle(bool flag){_useTriangle=flag;}
    void SetUseRayDiffEpsilon(bool flag){_useRayDiffEpsilon=flag;}
    void SetConservativeTest(bool flag){_conservativeTest=flag;}
    void SetDirectBilinear(bool flag){_directBilinear=flag;}

    const std::vector<BVHNode> &GetNodes() const { return _nodes; }
    const std::vector<unsigned int> &GetIndices() const { return _indices; }

private:
    ///< Builds BVH tree recursively.
    size_t BuildTree(const Mesh *mesh, unsigned int leftIdx,
                     unsigned int rightIdx, int depth);

    BVHBuildOptions _options;
    std::vector<BVHNode> _nodes;
    std::vector<unsigned int> _indices; // max 4G triangles.
    BVHBuildStatistics _stats;

    int _intersectKernel;
    float _uvMargin;
    bool _cropUV;
    bool _bezierClip;
    float _displaceScale;
    float _displaceFreq;

    double _epsilon;
    int    _maxLevel; 
    bool _useTriangle;
    bool _useRayDiffEpsilon;
    bool _conservativeTest;
    bool _directBilinear;
};

#endif // __BVH_ACCEL_H__
       // vim:set sw=2 ts=2 expandtab:
