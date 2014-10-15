#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include <vector>
#include <far/patchTables.h>
#include "common.h"

struct Mesh {
    Mesh() : numTriangles(0),
             numBezierPatches(0), displaceBound(0) {
    }

    bool IsBezierMesh() const {
        return numTriangles == 0;
    }

    // triangles
    size_t numTriangles;
    std::vector<real> triVertices;
    std::vector<real> triNormals;
    std::vector<unsigned int> faces;

    // patches
    size_t numBezierPatches;
    std::vector<real> bezierVertices;              /// [xyz] * 16 * numBezierPatches
    std::vector<real> bezierBounds;                /// [xyz] * [min, max] * numBezierPatches

    OpenSubdiv::Far::PatchTables::PatchParamTable patchParams;       /// [FarPatchParam] * numBezierPatches
    std::vector<float> colors;                     /// [rgb] * numBezierPatches;

    std::vector<int> wcpFlags; 


    float displaceBound;
};

#endif // __MESH_H__
