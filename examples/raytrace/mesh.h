#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include <vector>
#include <far/patchParam.h>
#include "common.h"

typedef struct {
    size_t numBezierPatches;
    std::vector<real> bezierVertices;              /// [xyz] * 16 * numBezierPatches
    OpenSubdiv::FarPatchParam const *patchParams;  /// [FarPatchParam] * numBezierPatches
    std::vector<float> colors;                     /// [rgb] * numBezierPatches;

} Mesh;

#endif // __MESH_H__
