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
                    OpenSubdiv::Far::PatchTables const *patchTables,
                    bool watertight,
                    float displaceBound)
{
    using namespace OpenSubdiv;

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
            vertices.push_back((v[0]-center[0])/radius);
            vertices.push_back((v[1]-center[1])/radius);
            //vertices.push_back(0);
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
    _wcpFlags.clear();

    std::vector<int> cpIndices; // 16 * numPatches

    int gregoryBegin = -1, gregoryEnd = 0;
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
            gregoryBegin = numTotalPatches;
            gregoryEnd = numTotalPatches + numPatches;
            break;
        case Far::PatchTables::GREGORY_BOUNDARY:
            numPatches = convertBoundaryGregory(_bezierVertices,
                                                _bezierBounds,
                                                cpIndices,
                                                &vertices[0], patchTables, *it);
            if (gregoryBegin == -1) gregoryBegin = numTotalPatches;
            gregoryEnd = numTotalPatches + numPatches;
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

        if (desc.GetType() == Far::PatchTables::SINGLE_CREASE) {
            for (int i = 0; i < numPatches; ++i) {
                patchParam.push_back(srcPatchParam[it->GetPatchIndex() + i]);
                patchParam.push_back(srcPatchParam[it->GetPatchIndex() + i]);
            }
        } else {
            for (int i = 0; i < numPatches; ++i) {
                patchParam.push_back(srcPatchParam[it->GetPatchIndex() + i]);
            }
        }
    }

    printf("%d %d\n", (int)patchParam.size(), (int)srcPatchParam.size());

    // for (int i = 0; i < hbrMesh->GetNumFaces(); ++i) {
    //     OsdHbrFace *f = hbrMesh->GetFace(i);
    //     bool out = ((f->isTransitionPatch() or f->_adaptiveFlags.patchType!=OsdHbrFace::kEnd) and
    //                 (not f->_adaptiveFlags.isExtraordinary) and
    //                 (f->_adaptiveFlags.bverts!=1));
    //     VERBOSE("Patch %d out=%x\n", i, out);
    // }

    _wcpFlags.resize(numTotalPatches);

    // vertex position verification pass
    if (watertight) {
        using namespace OsdBezier;

        std::map<Edge, Bezier > edgeBeziers;

        for (int i = 0; i < numTotalPatches; ++i) {
//            VERBOSE("\n============ patch %d (HBR : %d)==============\n", i, hbrMesh->GetFace(pathcParam[i].hbrFaceIndex)->GetID());
            int parentQuad[] = { 0, 2, 1, 3 };

            int edgeVerts[][4] = { {0, 1, 2, 3}, {0, 4, 8, 12}, {3, 7, 11, 15}, {12, 13, 14, 15} };
//            int edgeVerts[][4] = { {0, 4, 8, 12},{12, 13, 14, 15}, {15, 11, 7, 3}, {3, 2, 1, 0} };

            int edgeParents[][2] = { {0, 2}, {0, 1}, {2, 3}, {1, 3} };

            // store bezier edges (skip gregory)
            if (gregoryBegin <= i && i < gregoryEnd) continue;

            for (int j = 0; j < 4; ++j) {
                int parent = cpIndices[i*4+parentQuad[j]];
                if (parent < 0) continue;
                // parent mapping
                // parent = vertexParentIDs[parent];
                //NO_HBRparent = farToHbr[parent];
                //VERBOSE("%d/%d  %d: %f %f %f\n", i, cornerQuad[j], parent, v[0], v[1], v[2]);

                // for edge verts
                vec3f e0(&_bezierVertices[(i*16 + edgeVerts[j][0]) * 3]);
                vec3f e1(&_bezierVertices[(i*16 + edgeVerts[j][1]) * 3]);
                vec3f e2(&_bezierVertices[(i*16 + edgeVerts[j][2]) * 3]);
                vec3f e3(&_bezierVertices[(i*16 + edgeVerts[j][3]) * 3]);

                int ep0 = cpIndices[i*4 + edgeParents[j][0]];
                int ep1 = cpIndices[i*4 + edgeParents[j][1]];
                //ep0 = vertexParentIDs[ep0];
                //ep1 = vertexParentIDs[ep1];
//NO_HBR                ep0 = farToHbr[ep0];
//                ep1 = farToHbr[ep1];
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
        // save wcp flags

#if 0 // NO_HBR
        for (int i = 0; i < numTotalPatches; ++i) {
            int hbrFace = patchParam[i].hbrFaceIndex;
            OsdHbrFace *face = hbrMesh->GetFace(hbrFace);

            int wcpFlag = 0;
            int rots = face->_adaptiveFlags.rots;
            switch(face->_adaptiveFlags.transitionType) {
            case OsdHbrFace::kTransition0:
                wcpFlag = 1 << rots; break;
            case OsdHbrFace::kTransition1:
                wcpFlag = (1 << rots) | (1 << ((rots+3)%4)); break;
            case OsdHbrFace::kTransition2:
                wcpFlag = (1 << ((rots+1)%4)) | (1 << ((rots+2)%4)) | (1 << ((rots+3)%4));
            case OsdHbrFace::kTransition3:
                wcpFlag = 0xf;
            case OsdHbrFace::kTransition4:
                wcpFlag = (1 << rots) | (1 << ((rots+2)%4)); break;
                break;
            }
            wcpFlags[i] = wcpFlag;
        }
#endif

#if 0 //NO_HBR
        // different level subdivision
        // water-tight critical
        for (int i = 0; i < numTotalPatches; ++i) {
            int hbrFace = patchParam[i].hbrFaceIndex;
            OsdHbrFace *face = hbrMesh->GetFace(hbrFace);

            // ---------------- gregory
            if (gregoryBegin <= i && i < gregoryEnd) {
                for (int j = 0; j < 4; ++j) {
                    OsdHbrHalfedge *ed = face->GetEdge(j);
                    // pick two vertices
                    OsdHbrVertex *v0 = ed->GetVertex();
                    OsdHbrVertex *v1 = ed->GetNext()->GetVertex();

                    Edge e(v0->GetID(), v1->GetID());
                    // lookup edge dictionary
                    if (edgeBeziers.count(e) != 0) {
                        Bezier edge = edgeBeziers[e];
                        bool reverse = e._edges[0] != v0->GetID();

                        // overwrite
                        vec3f *d = (vec3f*)(&bezierVertices[i*16*3]);
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
                    } else {
                        // boundary of gregory-and-gregory.
                    }
                }
                continue;
            }

            // ---------------- non gregory

            if (not face->_adaptiveFlags.isCritical) continue;

            VERBOSE("Critical patch %d, %d----\n", i, hbrFace);
            OsdHbrFace *parentFace = face->GetParent();

            if (not parentFace) continue;

            if (parentFace->GetNumVertices() != 4) {
                VERBOSE("non-quad parent %d\n", parentFace->GetID());
                continue;
            }

            int childIndex = -1;
            bool watertightEdges[4] = { false, false, false, false };
            for (int j = 0; j < 4; ++j) {
                OsdHbrHalfedge *edge = parentFace->GetEdge(j)->GetOpposite();
                if (edge) {
                    OsdHbrFace *f = edge->GetFace();
                    
                    bool out = ((f->isTransitionPatch() or f->_adaptiveFlags.patchType!=OsdHbrFace::kEnd) and
                                (not f->_adaptiveFlags.isExtraordinary) and
                                (f->_adaptiveFlags.bverts!=1));

                    watertightEdges[j] = out;
                    VERBOSE("[%d: face=%d end=%d] ", j, f->GetID(), out);
                }
                if(parentFace->GetChild(j) == face) childIndex = j;
            }
            VERBOSE(" CHILD=%d, rot=%d\n", childIndex, face->_adaptiveFlags.rots);

            if (childIndex == -1) continue;

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

            int edgeScan[4][2] = { {3, 0}, {0, 1}, {1, 2}, {2, 3} };
            for (int j = 0; j < 2; ++j) {
                int edge = edgeScan[childIndex][j];
                if (not watertightEdges[edge]) continue;

                // at this point, we know the parent edge is T-node

                OsdHbrHalfedge *ed = parentFace->GetEdge(edge);

                // pick two vertices of the parent face
                OsdHbrVertex *v0 = ed->GetVertex();//parentFace->GetVertex(edge);
                OsdHbrVertex *v1 = ed->GetNext()->GetVertex(); //parentFace->GetVertex((edge+1)%4);

                Edge e(v0->GetID(), v1->GetID());
                VERBOSE( " search %d-%d edge, face %d \n", v0->GetID(), v1->GetID(),
                         parentFace->GetEdge(edge)->GetOpposite()->GetFace()->GetID());

                if (parentFace->_adaptiveFlags.bverts > 0) continue;

                if (edgeBeziers.count(e) != 0) {
                    Bezier parentEdge = edgeBeziers[e];  // todo, lookup
                    bool reverse = e._edges[0] != v0->GetID();
                    if (reverse) {
                        parentEdge.Reverse();
                    }

                    // find corresponding edge at patch[i]
                    // cut original bezier
                    vec3f tmp[2][4];
                    OsdBezier::BezierSplit<4, 0, /*stride*/1, vec3f, float> split(
                        tmp[0], tmp[1], edgeBeziers[e].cp, 0.5);
                    VERBOSE("[[%f %f %f - ", tmp[0][0][0], tmp[0][0][1], tmp[0][0][2]);
                    VERBOSE("%f %f %f - ", tmp[1][0][0], tmp[1][0][1], tmp[1][0][2]);
                    VERBOSE("%f %f %f]]\n", tmp[1][3][0], tmp[1][3][1], tmp[1][3][2]);

                    vec3f *d = (vec3f*)(&bezierVertices[i*16*3]);
                    if (edge == 0) {
                        for (int k = 0; k < 4; ++k) {
                            d[k*4+0] = reverse
                                ? tmp[childIndex==0?1:0][3-k]
                                : tmp[childIndex==0?0:1][k];
                        }
                    } else if (edge == 1) {
                        for (int k = 0; k < 4; ++k) {
                            d[12+k] = reverse
                                ? tmp[childIndex==1?1:0][3-k]
                                : tmp[childIndex==1?0:1][k];
                        }
                    } else if (edge == 2) {
                        for (int k = 0; k < 4; ++k) {
                            d[k*4+3] = reverse
                                ? tmp[childIndex==2?1:0][k]
                                : tmp[childIndex==2?0:1][3-k];
                        }
                    } else if (edge == 3) {
                        for (int k = 0; k < 4; ++k) {
                            d[k] = reverse
                                ? tmp[childIndex==3?1:0][k]
                                : tmp[childIndex==3?0:1][3-k];
                        }
                    }
                } else {
                    printf("Topology error --- Not found hbr verts in edge dict%d,%d\n",
                           v0->GetID(), v1->GetID());
                }
            }
        }
#endif
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
    _colors.clear();
    _numTriangles = 0;

    if (level == 0) {
        return;
    }

    float *cp = &_bezierVertices[0];
    int vindex = 0;
    int numTriangles = 0;

    std::vector<float> colors;
    for (size_t patchIndex = 0; patchIndex < _numBezierPatches; ++patchIndex) {

        const OpenSubdiv::Far::PatchParam &param = _patchParams[patchIndex];
        unsigned int bits = param.bitField.field;
        int patchLevel = (bits & 0xf);
        int div = 1 << std::max(0, level-patchLevel);

        float p[3], n[3];
        for (int u = 0; u <= div; ++u) {
            for (int v = 0; v <= div; ++v) {
                evalBezier(p, n, u/float(div), v/float(div), cp);
                _triVertices.push_back(p[0]);
                _triVertices.push_back(p[1]);
                _triVertices.push_back(p[2]);
                _triNormals.push_back(n[0]);
                _triNormals.push_back(n[1]);
                _triNormals.push_back(n[2]);

                if (u < div and v < div) {
                    _faces.push_back(vindex);
                    _faces.push_back(vindex+1);
                    _faces.push_back(vindex+div+1);

                    _faces.push_back(vindex+div+1);
                    _faces.push_back(vindex+1);
                    _faces.push_back(vindex+div+2);

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
