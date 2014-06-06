#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include <vector>
#include <far/patchParam.h>
#include "common.h"

struct Mesh {
    Mesh() : numTriangles(0), numBezierPatches(0), patchParams(NULL) {
    }

    bool IsBezierMesh() const {
        return numBezierPatches > 0;
    }

    size_t numTriangles;
    std::vector<real> triVertices;
    std::vector<unsigned int> faces;

    size_t numBezierPatches;
    std::vector<real> bezierVertices;              /// [xyz] * 16 * numBezierPatches
    std::vector<real> bezierBounds;                /// [xyz] * [min, max] * numBezierPatches
    OpenSubdiv::FarPatchParam const *patchParams;  /// [FarPatchParam] * numBezierPatches
    std::vector<float> colors;                     /// [rgb] * numBezierPatches;
};

#endif // __MESH_H__
