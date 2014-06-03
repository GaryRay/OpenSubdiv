#ifndef CONVERT_BEZIER_H
#define CONVERT_BEZIER_H

#include <vector>
#include <far/patchTables.h>

extern int convertRegular(std::vector<float> &bezierVertices,
                          float const *vertices,
                          OpenSubdiv::FarPatchTables const *patchTables,
                          OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertBoundary(std::vector<float> &bezierVertices,
                          float const *vertices,
                           OpenSubdiv::FarPatchTables const *patchTables,
                           OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertCorner(std::vector<float> &bezierVertices,
                          float const *vertices,
                         OpenSubdiv::FarPatchTables const *patchTables,
                         OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertGregory(std::vector<float> &bezierVertices,
                          float const *vertices,
                          OpenSubdiv::FarPatchTables const *patchTables,
                          OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertBoundaryGregory(std::vector<float> &bezierVertices,
                                  float const *vertices,
                                  OpenSubdiv::FarPatchTables const *patchTables,
                                  OpenSubdiv::FarPatchTables::PatchArray const &parray);


#endif  // CONVERT_BEZIER_H
