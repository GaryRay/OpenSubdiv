#ifndef CONVERT_BEZIER_H
#define CONVERT_BEZIER_H

#include <vector>
#include <far/patchTables.h>

extern int convertRegular(std::vector<float> &bezierVertices,
                          std::vector<float> &bezierBounds,
                          std::vector<int> &cpIndices,
                          float const *vertices,
                          OpenSubdiv::Far::PatchTables const *patchTables,
                          OpenSubdiv::Far::PatchTables::PatchArray const &parray);

extern int convertSingleCrease(std::vector<float> &bezierVertices,
                               std::vector<float> &bezierBounds,
                               std::vector<int> &cpIndices,
                               float const *vertices,
                               OpenSubdiv::Far::PatchTables const *patchTables,
                               OpenSubdiv::Far::PatchTables::PatchArray const &parray);

extern int convertBoundary(std::vector<float> &bezierVertices,
                           std::vector<float> &bezierBounds,
                           std::vector<int> &cpIndices,
                           float const *vertices,
                           OpenSubdiv::Far::PatchTables const *patchTables,
                           OpenSubdiv::Far::PatchTables::PatchArray const &parray);

extern int convertCorner(std::vector<float> &bezierVertices,
                         std::vector<float> &bezierBounds,
                         std::vector<int> &cpIndices,
                         float const *vertices,
                         OpenSubdiv::Far::PatchTables const *patchTables,
                         OpenSubdiv::Far::PatchTables::PatchArray const &parray);

extern int convertGregory(std::vector<float> &bezierVertices,
                          std::vector<float> &bezierBounds,
                          std::vector<int> &cpIndices,
                          float const *vertices,
                          OpenSubdiv::Far::PatchTables const *patchTables,
                          OpenSubdiv::Far::PatchTables::PatchArray const &parray);

extern int convertBoundaryGregory(std::vector<float> &bezierVertices,
                                  std::vector<float> &bezierBounds,
                                  std::vector<int> &cpIndices,
                                  float const *vertices,
                                  OpenSubdiv::Far::PatchTables const *patchTables,
                                  OpenSubdiv::Far::PatchTables::PatchArray const &parray);


#endif  // CONVERT_BEZIER_H
