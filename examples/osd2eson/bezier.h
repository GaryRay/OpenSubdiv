#ifndef BEZIER_H
#define BEZIER_H

#include <vector>
#include <far/patchTables.h>

extern int convertRegular(std::vector<float> &bezierVertices,
                          std::vector<float> const &vertices,
                          OpenSubdiv::FarPatchTables const *patchTables,
                          OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertBoundary(std::vector<float> &bezierVertices,
                           std::vector<float> const &vertices,
                           OpenSubdiv::FarPatchTables const *patchTables,
                           OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertCorner(std::vector<float> &bezierVertices,
                         std::vector<float> const &vertices,
                         OpenSubdiv::FarPatchTables const *patchTables,
                         OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertGregory(std::vector<float> &bezierVertices,
                          std::vector<float> const &vertices,
                          OpenSubdiv::FarPatchTables const *patchTables,
                          OpenSubdiv::FarPatchTables::PatchArray const &parray);

extern int convertBoundaryGregory(std::vector<float> &bezierVertices,
                                  std::vector<float> const &vertices,
                                  OpenSubdiv::FarPatchTables const *patchTables,
                                  OpenSubdiv::FarPatchTables::PatchArray const &parray);


#endif  // BEZIER_H
