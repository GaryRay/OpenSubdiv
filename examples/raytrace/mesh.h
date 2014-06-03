#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include "common.h"

#include <far/patchParam.h>

typedef struct {
  size_t numVertices;
  real *vertices;              /// [xyz] * numVertices
  size_t numBezierPatches;
  real *bezierVertices;        /// [xyz] * 16 * numBezierPatches
  OpenSubdiv::FarPatchParam const *patchParams;  /// [FarPatchParam] * numBezierPatches

} Mesh;

#endif // __MESH_H__
