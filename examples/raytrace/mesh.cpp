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

#include <cfloat>
#include <set>
#include "mesh.h"
#include "convert_bezier.h"
#include "bezier/math.h"
#include "bezier/bezier.h"
#include "../common/patchColors.h"

#define VERBOSE(x, ...)
//#define VERBOSE(x, ...) printf(x, __VA_ARGS__)

struct Edge {
    Edge(int a, int b) {
        if (a > b) std::swap(a, b);
        _edges[0] = a;
        _edges[1] = b;
    }
    bool operator < (Edge const &other) const {
        return (_edges[0] < other._edges[0] ||
                (_edges[0] == other._edges[0] && _edges[1]  < other._edges[1]));
    }
    int _edges[2];
};

struct Bezier {
    Bezier() {}
    Bezier(OsdBezier::vec3f p[4]) {
        cp[0] = p[0];
        cp[1] = p[1];
        cp[2] = p[2];
        cp[3] = p[3];
    }
    Bezier(OsdBezier::vec3f const &p0, OsdBezier::vec3f const &p1,
           OsdBezier::vec3f const &p2, OsdBezier::vec3f const &p3) {
        cp[0] = p0;
        cp[1] = p1;
        cp[2] = p2;
        cp[3] = p3;
    }
    void Reverse() {
        std::swap(cp[0], cp[3]);
        std::swap(cp[1], cp[2]);
    }
    OsdBezier::vec3f cp[4];
};

void
Mesh::BezierConvert(float *inVertices, int numVertices,
                    OpenSubdiv::Far::TopologyRefiner const *refiner,
                    OpenSubdiv::Far::PatchTables const *patchTables,
                    bool watertight,
                    float displaceBound)
{
    using namespace OpenSubdiv;
    using namespace OsdBezier;

    // convert to mesh
    Far::PatchTables::PatchArrayVector const &patchArrays =
        patchTables->GetPatchArrayVector();
    Far::PatchTables::PatchParamTable const &srcPatchParam =
        patchTables->GetPatchParamTable();

    Far::PatchTables::PatchParamTable patchParam;

    // centering/normalize vertices.
    std::vector<float> vertices;
    vertices.reserve(numVertices*3);
    {
        float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
        float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (int i = 0; i < numVertices; ++i) {
            float *v = inVertices + i*3;
            for (int j = 0; j < 3; ++j) {
                min[j] = std::min(min[j], v[j]);
                max[j] = std::max(max[j], v[j]);
            }
        }
        float center[3] = { (max[0]+min[0])*0.5f,
                            (max[1]+min[1])*0.5f,
                            (max[2]+min[2])*0.5f };
        float radius = std::max(std::max(max[0]-min[0], max[1]-min[1]), max[2]-min[2]);
        for (int i = 0; i < numVertices; ++i) {
            float *v = inVertices + i*3;
#if 1
            // centering
            vertices.push_back((v[0]-center[0])/radius);
            vertices.push_back((v[1]-center[1])/radius);
            vertices.push_back((v[2]-center[2])/radius);
#else
            vertices.push_back(v[0]);
            vertices.push_back(v[1]);
            vertices.push_back(v[2]);
#endif
        }
    }

    int numTotalPatches = 0;
    _bezierVertices.clear();
    _bezierBounds.clear();
    _colors.clear();
    _sharpnesses.clear();
    _wcpFlags.clear();

    std::vector<int> cpIndices; // 16 * numPatches

    std::vector<Far::PatchTables::Descriptor> patchDescs;

    // iterate patch types.
    for (Far::PatchTables::PatchArrayVector::const_iterator it = patchArrays.begin();
         it != patchArrays.end(); ++it) {

        int numPatches = 0;
        Far::PatchTables::Descriptor desc = it->GetDescriptor();
        VERBOSE("TransitionType = %d\n", desc.GetPattern());

        switch(desc.GetType()) {
        case Far::PatchTables::REGULAR:
            numPatches = convertRegular(_bezierVertices,
                                        _bezierBounds,
                                        cpIndices,
                                        &vertices[0], patchTables, *it);
            break;
        case Far::PatchTables::SINGLE_CREASE:
            numPatches = convertSingleCrease(_bezierVertices,
                                             _bezierBounds,
                                             cpIndices,
                                             &vertices[0], patchTables, *it);
            break;
        case Far::PatchTables::BOUNDARY:
            numPatches = convertBoundary(_bezierVertices,
                                         _bezierBounds,
                                         cpIndices,
                                         &vertices[0], patchTables, *it);
            break;
        case Far::PatchTables::CORNER:
            numPatches = convertCorner(_bezierVertices,
                                       _bezierBounds,
                                       cpIndices,
                                       &vertices[0], patchTables, *it);
            break;
        case Far::PatchTables::GREGORY:
            numPatches = convertGregory(_bezierVertices,
                                        _bezierBounds,
                                        cpIndices,
                                        &vertices[0], patchTables, *it);
            break;
        case Far::PatchTables::GREGORY_BOUNDARY:
            numPatches = convertBoundaryGregory(_bezierVertices,
                                                _bezierBounds,
                                                cpIndices,
                                                &vertices[0], patchTables, *it);
            break;
        default:
            break;
        }
        numTotalPatches += numPatches;
        // color array
        {
            const float *color = getAdaptivePatchColor(desc);
            for (int i = 0; i < numPatches; ++i) {
                _colors.push_back(color[0]);
                _colors.push_back(color[1]);
                _colors.push_back(color[2]);
            }
        }
        //descs
        {
            for (int i = 0; i < numPatches; ++i) {
                patchDescs.push_back(desc);
            }
        }

        if (desc.GetType() == Far::PatchTables::SINGLE_CREASE) {
            // duplicate patchparam
            for (int i = 0; i < numPatches/2; ++i) {
                patchParam.push_back(srcPatchParam[it->GetPatchIndex() + i]);
                patchParam.push_back(srcPatchParam[it->GetPatchIndex() + i]);
                int sharpid = patchTables->GetSharpnessIndexTable()[it->GetPatchIndex() + i];
                float sharpness = patchTables->GetSharpnessValues()[sharpid];
                _sharpnesses.push_back(sharpness);
                _sharpnesses.push_back(-sharpness);
            }
        } else {
            for (int i = 0; i < numPatches; ++i) {
                patchParam.push_back(srcPatchParam[it->GetPatchIndex() + i]);
                _sharpnesses.push_back(0);
            }
        }
    }

    _wcpFlags.resize(numTotalPatches);

    // vertex position verification pass
    if (watertight) {
        std::map<Edge, Bezier> edgeBeziers;

        // save wcp flags
        for (int i = 0; i < numTotalPatches; ++i) {
            int face = patchParam[i].faceIndex;
            int wcpFlag = 0;
            int rots = 0;//patchDescs[i].rotation??
            //int rots = face->_adaptiveFlags.rots;

            switch(patchDescs[i].GetPattern()) {
            case Far::PatchTables::PATTERN0:
                wcpFlag = 1 << rots; break;
                break;
            case Far::PatchTables::PATTERN1:
                wcpFlag = (1 << rots) | (1 << ((rots+3)%4)); break;
                break;
            case Far::PatchTables::PATTERN2:
                wcpFlag = (1 << ((rots+1)%4)) | (1 << ((rots+2)%4)) | (1 << ((rots+3)%4));
                break;
            case Far::PatchTables::PATTERN3:
                wcpFlag = 0xf;
                break;
            case Far::PatchTables::PATTERN4:
                wcpFlag = (1 << rots) | (1 << ((rots+2)%4)); break;
                break;
            }
            _wcpFlags[i] = wcpFlag;
        }


        for (int i = 0; i < numTotalPatches; ++i) {
            int parentQuad[] = { 0, 2, 1, 3 };

            int edgeVerts[][4] = { {0, 1, 2, 3}, {0, 4, 8, 12}, {3, 7, 11, 15}, {12, 13, 14, 15} };
            int edgeParents[][2] = { {0, 2}, {0, 1}, {2, 3}, {1, 3} };

            // store bezier edges (skip gregory, single-crease)
            if (patchDescs[i].GetType() == Far::PatchTables::SINGLE_CREASE or
                patchDescs[i].GetType() == Far::PatchTables::GREGORY or
                patchDescs[i].GetType() == Far::PatchTables::GREGORY_BOUNDARY) continue;

            VERBOSE("\n============ patch %d ==============\n", patchParam[i].faceIndex);

            for (int j = 0; j < 4; ++j) {
                int parent = cpIndices[i*4+parentQuad[j]];
                if (parent < 0) continue; // skip border, boundary (for now)
                // parent mapping
                // parent = vertexParentIDs[parent];
                //NO_HBRparent = farToHbr[parent];
                //VERBOSE("%d  %d: %f %f %f\n", i, parent, v[0], v[1], v[2]);

                // for edge verts
                vec3f e0(&_bezierVertices[(i*16 + edgeVerts[j][0]) * 3]);
                vec3f e1(&_bezierVertices[(i*16 + edgeVerts[j][1]) * 3]);
                vec3f e2(&_bezierVertices[(i*16 + edgeVerts[j][2]) * 3]);
                vec3f e3(&_bezierVertices[(i*16 + edgeVerts[j][3]) * 3]);

                int ep0 = cpIndices[i*4 + edgeParents[j][0]];
                int ep1 = cpIndices[i*4 + edgeParents[j][1]];
                Edge edge(ep0, ep1);
                VERBOSE("(%d, %d) (%f, %f, %f)-(%f, %f, %f)\n", ep0, ep1,
                       e0[0], e0[1], e0[2],
                       e3[0], e3[1], e3[2]);
                if (edgeBeziers.count(edge) == 0) {
                    if (edge._edges[0] == ep0) {
                        edgeBeziers[edge] = Bezier(e0, e1, e2, e3);
                    } else {
                        edgeBeziers[edge] = Bezier(e3, e2, e1, e0);
                    }
                }
            }
        }
  /*
    pattern0        pattern1       pattern2        pattern3        pattern4
 +-------------+ +-------------+ +-------------+ +-------------+ +-------------+
 |     /\\     | | 0   /\\   2 | |             | |      |      | |      |      |
 | 1  /  \\  2 | |    /   \\   | |      0      | |  1   |  0   | |      |      |
 |   /    \\   | |   /  3   \\ | |-------------| |------|------| |  1   |   0  |
 |  /      \\  | |  /       /  | |\\    3    / | |      |      | |      |      |
 | /    0   \\ | | /    /    1 | |  \\     /   | |  3   |  2   | |      |      |
 |/          \\| |/ /          | | 1  \\ /   2 | |      |      | |      |      |
 +-------------+ +-------------+ +-------------+ +-------------+ +-------------+
*/
        std::map<std::pair<int, int>, int> finalFaces;
        for (int i = 0; i < numTotalPatches; ++i) {
            int faceIndex = patchParam[i].vtrFaceIndex;
            int level = patchParam[i].bitField.GetDepth();
            finalFaces[std::make_pair(level, faceIndex)] = i; // what if single-crease?
        }

        // different level subdivision
        // water-tight critical
        for (int i = 0; i < numTotalPatches; ++i) {
            //int hbrFace = patchParam[i].hbrFaceIndex;
            //OsdHbrFace *face = hbrMesh->GetFace(hbrFace);
            int faceIndex = patchParam[i].vtrFaceIndex;

            if (patchDescs[i].GetType() == Far::PatchTables::SINGLE_CREASE) continue;
            // ---------------- gregory
            //   0-----3
            //   |     |
            //   |     |
            //   1-----2
            if (patchDescs[i].GetType() == Far::PatchTables::GREGORY or
                patchDescs[i].GetType() == Far::PatchTables::GREGORY_BOUNDARY) {
                int edgeParents[][2] = { {0, 1}, {1, 2}, {2, 3}, {3, 0} };
                for (int j = 0; j < 4; ++j) {
                    int v0 = cpIndices[i*4 + edgeParents[j][0]];
                    int v1 = cpIndices[i*4 + edgeParents[j][1]];
                    Edge e(v0, v1);
                    if (edgeBeziers.count(e) != 0) {
                        // overwrite
                        Bezier edge = edgeBeziers[e];
                        bool reverse = e._edges[0] != v0;
                        vec3f *d = (vec3f*)(&_bezierVertices[i*16*3]);
                        if (j == 0) {
                            for (int k = 0; k < 4; ++k) {
                                d[k*4+0] = reverse ? edge.cp[3-k] : edge.cp[k];
                            }
                        } else if (j == 1) {
                            for (int k = 0; k < 4; ++k) {
                                d[12+k] = reverse ? edge.cp[3-k] : edge.cp[k];
                            }
                        } else if (j == 2) {
                            for (int k = 0; k < 4; ++k) {
                                d[k*4+3] = reverse ? edge.cp[k] : edge.cp[3-k];
                            }
                        } else if (j == 3) {
                            for (int k = 0; k < 4; ++k) {
                                d[k] = reverse ? edge.cp[k] : edge.cp[3-k];
                            }
                        }
                    }
                }
                continue;
            }

            // ---------------- non gregory
            //if (_wcpFlags[i] == 0) continue;
            VERBOSE("\n");

            int level = patchParam[i].bitField.GetDepth()-1;
            if (level < 0) continue;
            int parentFace = refiner->GetChildFaceParentFace(level, faceIndex);
            Vtr::IndexArray indices = refiner->GetFaceVertices(level, parentFace);

            VERBOSE("Level %d : \e[32mPatch %d\e[0m, parent = %d (Flag=%d)----\n",
                    level+1, faceIndex, parentFace, _wcpFlags[i]);

            Vtr::IndexArray myEdges = refiner->GetFaceEdges(level+1, faceIndex);
            VERBOSE("My edges verts:");
            for (int i = 0; i < myEdges.size(); ++i) {
                VERBOSE(" %d ", myEdges[i]);
            }
            VERBOSE("\n");

            VERBOSE("Parent verts: level(%d)= ", level);
            for (int i = 0; i < indices.size(); ++i) {
                VERBOSE(" %d ", indices[i]);
            }
            VERBOSE("\n");

            if (indices.size() != 4) {
                VERBOSE("non-quad parent (nverts=%d)\n", (int)indices.size());
                continue;
            }

            // locate which child this face is, in the parent face (within 0-3)

            int childIndex = -1;
            bool watertightEdges[4] = { false, false, false, false };
            Vtr::IndexArray children = refiner->GetFaceChildFaces(level, parentFace);
            Vtr::IndexArray edges = refiner->GetFaceEdges(level, parentFace);

            for (int j = 0; j < 4; ++j) {
                if (children[j] == faceIndex) childIndex = j;
            }
            //int rot = patchDescs[i].GetRotation();
            int rot = patchParam[i].bitField.GetRotation();
            VERBOSE(" CHILD=%d, \e[31mrot=%d\e[0m\n", childIndex, rot);//face->_adaptiveFlags.rots);
            if (childIndex == -1) continue;

            VERBOSE("Child faces: (childIndex = %d)", childIndex);
            for (int i = 0; i < children.size(); ++i) {
                VERBOSE(" %d ", children[i]);
            }
            VERBOSE("\n");

            // check 2 edges.
            for (int j = 0; j < 2; ++j) {
                int edgeIndex = (j+3+childIndex)%4;

                //edgeIndex = (edgeIndex + rot)%4;

                if (patchDescs[i].GetType() == Far::PatchTables::BOUNDARY or
                    patchDescs[i].GetType() == Far::PatchTables::CORNER) continue;

                // if it's triangle head, skip
                //if ((_wcpFlags[i] >> edgeIndex)&1 == 1) continue;
                //if (_wcpFlags[i]) continue;
                
                Vtr::IndexArray faces = refiner->GetEdgeFaces(level, edges[edgeIndex]);
                if (faces.size() != 2) continue;
                int adjFace = (faces[0] == parentFace) ? faces[1] : faces[0];
                VERBOSE("Adjacent face %d\n", adjFace);

                // see if the adjacent face is final or not
                if (finalFaces.count(std::make_pair(level, adjFace)) > 0) {
                    VERBOSE("Critical Patch (level=%d, adjFace=%d exists)!\n", level, adjFace);
                    watertightEdges[edgeIndex] = true;
                } else {
                    continue;
                }

                int patchIndex = finalFaces[std::make_pair(level, adjFace)];


                // inefficient!
                int levelVertsOffset = 0;
                for (int k = 0; k < level; ++k){
                    levelVertsOffset += refiner->GetNumVertices(k);
                }
                int v0 = indices[edgeIndex] + levelVertsOffset;
                int v1 = indices[(edgeIndex+1)%4] + levelVertsOffset;
                VERBOSE("Matching verts: (level %d), %d- %d,  edge = %d\n", level,
                       v0, v1, edges[edgeIndex]);

                Edge e(v0, v1);

                //if (parentFace->_adaptiveFlags.bverts > 0) continue;
                // childindex 0 -> edge 3, 0
                // childindex 1 -> edge 0, 1
                // childindex 2 -> edge 1, 2
                // childindex 3 -> edge 2, 3
                // for each edge

                    /*
                      <----------- v
                          edge 0
                      +-----+-----+         u
                      |     |     |         |
                      |  1  |  0  |         |
                      |     |     |         |
               edge 1 +-----+-----+ edge 3  |
                      |     |     |         |
                      |  2  |  3  |         v
                      |     |(childIndex)|
                      +-----+-----+
                          edge 2
                     */

                if (edgeBeziers.count(e) != 0) {
                    Bezier parentEdge = edgeBeziers[e];  // todo, fix lookup
                    bool reverse = e._edges[0] != v0;
                    if (reverse) {
                        parentEdge.Reverse();
                    }

                    // apply rotation (for transition patches)
                    int edge = edgeIndex;
                    edge = (edge+(4-rot))%4;

                    VERBOSE("Align edgeIndex = %d, edge = %d\n", edgeIndex, edge);

                    // find corresponding edge at patch[i]
                    // cut original bezier
                    vec3f tmp[2][4];
                    OsdBezier::BezierSplit<4, 0, /*stride*/1, vec3f, float> split(
                        tmp[0], tmp[1], edgeBeziers[e].cp, 0.5);
                    VERBOSE("[[%f %f %f - ", tmp[0][0][0], tmp[0][0][1], tmp[0][0][2]);
                    VERBOSE("%f %f %f - ", tmp[1][0][0], tmp[1][0][1], tmp[1][0][2]);
                    VERBOSE("%f %f %f]]\n", tmp[1][3][0], tmp[1][3][1], tmp[1][3][2]);

                    vec3f *d = (vec3f*)(&_bezierVertices[i*16*3]);
                    if (edge == 0) {
                        for (int k = 0; k < 4; ++k) {
                            d[k*4+0] = reverse
                                ? tmp[childIndex==edgeIndex?1:0][3-k]
                                : tmp[childIndex==edgeIndex?0:1][k];
                        }
                    } else if (edge == 1) {
                        for (int k = 0; k < 4; ++k) {
                            d[12+k] = reverse
                                ? tmp[childIndex==edgeIndex?1:0][3-k]
                                : tmp[childIndex==edgeIndex?0:1][k];
                        }
                    } else if (edge == 2) {
                        for (int k = 0; k < 4; ++k) {
                            d[k*4+3] = reverse
                                ? tmp[childIndex==edgeIndex?1:0][k]
                                : tmp[childIndex==edgeIndex?0:1][3-k];
                        }
                    } else if (edge == 3) {
                        for (int k = 0; k < 4; ++k) {
                            d[k] = reverse
                                ? tmp[childIndex==edgeIndex?1:0][k]
                                : tmp[childIndex==edgeIndex?0:1][3-k];
                        }
                    }
                } else {
                    // not found in the edge dictionary.
                }
            }
        }
    }

    _numTriangles = 0;
    _numBezierPatches = numTotalPatches;
    _patchParams = patchParam;

    _displaceBound = displaceBound;

    assert(numTotalPatches*16*3 == (int)_bezierVertices.size());
}

static void evalBezier(float *p, float *n, float u, float v, const float *cp)
{
    OsdBezier::BezierPatch<OsdBezier::vec3f, float, 4> patch((const OsdBezier::vec3f*)cp);
    OsdBezier::vec3f b = patch.Evaluate(u, v);
    p[0] = b[0];
    p[1] = b[1];
    p[2] = b[2];
    OsdBezier::vec3f du = patch.EvaluateDu(u, v);
    OsdBezier::vec3f dv = patch.EvaluateDv(u, v);
    n[0] = -(du[1] * dv[2] - du[2] * dv[1]);
    n[1] = -(du[2] * dv[0] - du[0] * dv[2]);
    n[2] = -(du[0] * dv[1] - du[1] * dv[0]);
}

void
Mesh::Tessellate(int level)
{
    _triVertices.clear();
    _triNormals.clear();
    _faces.clear();
    //_colors.clear();
    _numTriangles = 0;

    // _triVertices.reserve(120000000*3*3);
    // _triNormals.reserve(120000000*3*3);
    // _faces.reserve(120000000*3*2);

    std::vector<float> colors;
    //colors.reserve(120000000*3*2);

    if (level == 0) {
        return;
    }

    float *cp = &_bezierVertices[0];
    int vindex = 0;
    int numTriangles = 0;

    for (size_t patchIndex = 0; patchIndex < _numBezierPatches; ++patchIndex) {

        const OpenSubdiv::Far::PatchParam &param = _patchParams[patchIndex];
        unsigned int bits = param.bitField.field;
        int patchLevel = (bits & 0xf);
        int div = 1 << std::max(0, level-patchLevel);
        int udiv = div;
        int vdiv = div;

        float sharp = _sharpnesses[patchIndex];
        if (sharp > 0) {
            udiv = 1 << std::max(0, level-patchLevel-int(sharp));
        } else if (sharp < 0) {
            udiv = div - (1 << std::max(0, level-patchLevel+int(sharp)));
            udiv = std::max(1, udiv);
        }

        float p[3], n[3];
        for (int u = 0; u <= udiv; ++u) {
            for (int v = 0; v <= vdiv; ++v) {
                evalBezier(p, n, u/float(udiv), v/float(vdiv), cp);
                _triVertices.push_back(p[0]);
                _triVertices.push_back(p[1]);
                _triVertices.push_back(p[2]);
                _triNormals.push_back(n[0]);
                _triNormals.push_back(n[1]);
                _triNormals.push_back(n[2]);

                if (u < udiv and v < vdiv) {
                    _faces.push_back(vindex);
                    _faces.push_back(vindex+1);
                    _faces.push_back(vindex+vdiv+1);

                    _faces.push_back(vindex+vdiv+1);
                    _faces.push_back(vindex+1);
                    _faces.push_back(vindex+vdiv+2);

                    colors.push_back(_colors[patchIndex*3+0]);
                    colors.push_back(_colors[patchIndex*3+1]);
                    colors.push_back(_colors[patchIndex*3+2]);
                    colors.push_back(_colors[patchIndex*3+0]);
                    colors.push_back(_colors[patchIndex*3+1]);
                    colors.push_back(_colors[patchIndex*3+2]);
                    numTriangles += 2;
                }
                ++vindex;
            }
        }

        // next
        cp += 16*3;
    }

    _numTriangles = numTriangles;
    _colors.swap(colors);
}
