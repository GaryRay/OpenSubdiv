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
#ifndef __MESH_H__
#define __MESH_H__

#include <cstdio>
#include <vector>
#include <far/topologyRefiner.h>
#include <far/patchTables.h>

struct Material {
    float ka[3];         // ambient
    float kd[3];         // diffuse
    float ks[3];         // specular
    float ns;            // specular exponent
    float ni;            // optical density (1.0=no refraction, glass=1.5)
    float sharpness;     // reflection sharpness
    float tf[3];         // transmission filter
    float d;             // dissolve factor (1.0 = opaque)
};

struct Mesh {
    Mesh() : _numTriangles(0),
             _numBezierPatches(0),
             _displaceBound(0) {
    }

    void BezierConvert(float *vertices, int numVertices,
                       OpenSubdiv::Far::TopologyRefiner const *refiner,
                       OpenSubdiv::Far::PatchTables const *patchTables,
                       bool watertight,
                       float displaceBound);

    void Tessellate(int level);

    bool IsBezierMesh() const { return _numTriangles == 0; }
    int GetNumPatches() const { return _numBezierPatches; }
    int GetNumTriangles() const { return _numTriangles; }

    size_t GetMemoryUsage() const {
        size_t mem = 0;
        if (IsBezierMesh()) {
            mem += _bezierVertices.size() * sizeof(float);  // cp
            mem += _numBezierPatches/2; // (4bit per patch);
            //mem += _mesh.bezierBounds.size() * sizeof(float);   // bounds
        } else {
            mem += _triVertices.size() * sizeof(float); // verts
            mem += _faces.size() * sizeof(unsigned int); // indices
        }
        return mem;
    }

    // triangles
    size_t _numTriangles;
    std::vector<float> _triVertices;
    std::vector<float> _triNormals;
    std::vector<unsigned int> _faces;

    // patches
    size_t _numBezierPatches;
    std::vector<float> _bezierVertices;              /// [xyz] * 16 * numBezierPatches
    std::vector<float> _bezierBounds;                /// [xyz] * [min, max] * numBezierPatches

    OpenSubdiv::Far::PatchParamTable _patchParams;       /// [FarPatchParam] * numBezierPatches
    std::vector<float> _colors;                     /// [rgb] * numBezierPatches;
    std::vector<float> _sharpnesses;                /// sharpness * numBezierPatches (for single-crease patch)

    std::vector<int> _wcpFlags; 
    std::vector<int> _materialIDs;
    std::vector<Material> _materials;

    float _displaceBound;
};

#endif // __MESH_H__
