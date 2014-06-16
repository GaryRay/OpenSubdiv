#define NEED_HBR_FACE_INDEX
#include "scene.h"
#include "convert_bezier.h"
#include "camera.h"

#include <ctime>
#include <cstring>
#include <string>
#include <cfloat>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

#ifdef OPENSUBDIV_HAS_TBB
#include <tbb/parallel_for.h>
#include <tbb/task_scheduler_init.h>
#endif

#include <osd/cpuEvalLimitContext.h>
#include <osd/cpuEvalLimitController.h>
#include "osdutil/bezier.h"
#include "osdutil/bezierIntersect.h"
#include "osdutil/math.h"

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
    Bezier(OsdUtil::vec3f p[4]) {
        cp[0] = p[0];
        cp[1] = p[1];
        cp[2] = p[2];
        cp[3] = p[3];
    }
    Bezier(OsdUtil::vec3f const &p0, OsdUtil::vec3f const &p1,
           OsdUtil::vec3f const &p2, OsdUtil::vec3f const &p3) {
        cp[0] = p0;
        cp[1] = p1;
        cp[2] = p2;
        cp[3] = p3;
    }
    void Reverse() {
        std::swap(cp[0], cp[3]);
        std::swap(cp[1], cp[2]);
    }
    OsdUtil::vec3f cp[4];
};

static float diffBezier(const Bezier &a, const Bezier &b)
{
    return (a.cp[0]-b.cp[0]).length()
        + (a.cp[1]-b.cp[1]).length()
        + (a.cp[2]-b.cp[2]).length()
        + (a.cp[3]-b.cp[3]).length();
}

// brute-force matching
static bool consolidateBezier(const Bezier &r0, const Bezier &r1, float *v)
{
    using namespace OsdUtil;

    int edgeVerts[][4] = { {0, 4, 8, 12},
                           {12, 13, 14, 15},
                           {15, 11, 7, 3},
                           {3, 2, 1, 0} };

    Bezier testBezier[4];
    testBezier[0] = r0;
    testBezier[1] = r1;
    testBezier[2] = r0;
    testBezier[3] = r1;
    testBezier[2].Reverse();
    testBezier[3].Reverse();

    float dmin = 1.0f;
    int imin = -1, jmin = -1;
    for (int i = 0; i < 4; ++i) {
        Bezier bezier(vec3f(v + edgeVerts[i][0]*3),
                      vec3f(v + edgeVerts[i][1]*3),
                      vec3f(v + edgeVerts[i][2]*3),
                      vec3f(v + edgeVerts[i][3]*3));

        for (int j = 0; j < 4; ++j) {
            float d = diffBezier(testBezier[j], bezier);
            if (d < dmin) {
                dmin = d;
                imin = i;
                jmin = j;
            }
        }
    }

    if (dmin < 0.00001f && imin >=0 && jmin >= 0) {
//        VERBOSE("DMIN = %.10f\n", dmin);
        for (int k = 0; k < 4; ++k) {
            for (int l = 0; l < 3; ++l) {
                v[edgeVerts[imin][k]*3+l] = testBezier[jmin].cp[k][l];
            }
        }
        return true;
    } else {
        return false;
    }
}

static float const *getAdaptivePatchColor(OpenSubdiv::FarPatchTables::Descriptor const & desc)
{
    static float _colors[4][5][4] = {{{1.0f,  1.0f,  1.0f,  1.0f},   // regular
                                      {0.8f,  0.0f,  0.0f,  1.0f},   // boundary
                                      {0.0f,  1.0f,  0.0f,  1.0f},   // corner
                                      {1.0f,  1.0f,  0.0f,  1.0f},   // gregory
                                      {1.0f,  0.5f,  0.0f,  1.0f}},  // gregory boundary

                                     {{0.0f,  1.0f,  1.0f,  1.0f},   // regular pattern 0
                                      {0.0f,  0.5f,  1.0f,  1.0f},   // regular pattern 1
                                      {0.0f,  0.5f,  0.5f,  1.0f},   // regular pattern 2
                                      {0.5f,  0.0f,  1.0f,  1.0f},   // regular pattern 3
                                      {1.0f,  0.5f,  1.0f,  1.0f}},  // regular pattern 4
 
                                     {{0.0f,  0.0f,  0.75f, 1.0f},   // boundary pattern 0
                                      {0.0f,  0.2f,  0.75f, 1.0f},   // boundary pattern 1
                                      {0.0f,  0.4f,  0.75f, 1.0f},   // boundary pattern 2
                                      {0.0f,  0.6f,  0.75f, 1.0f},   // boundary pattern 3
                                      {0.0f,  0.8f,  0.75f, 1.0f}},  // boundary pattern 4
 
                                     {{0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 0
                                      {0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 1
                                      {0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 2
                                      {0.25f, 0.25f, 0.25f, 1.0f},   // corner pattern 3
                                      {0.25f, 0.25f, 0.25f, 1.0f}}}; // corner pattern 4

    typedef OpenSubdiv::FarPatchTables FPT;

    if (desc.GetPattern()==FPT::NON_TRANSITION) {
        return _colors[0][(int)(desc.GetType()-FPT::REGULAR)];
    } else {
        return _colors[(int)(desc.GetType()-FPT::REGULAR)+1][(int)desc.GetPattern()-1];
    }
}

Scene::Scene() : _vbo(0), _watertight(false)
{
#ifdef OPENSUBDIV_HAS_TBB
    static tbb::task_scheduler_init init;
#endif
}

Scene::~Scene()
{
    if (_vbo) glDeleteBuffers(1, &_vbo);
}

void
Scene::BezierConvert(float *inVertices, int numVertices,
                     OpenSubdiv::FarPatchTables const *patchTables,
                     std::vector<int> const &farToHbr,
                     OsdHbrMesh *hbrMesh)
{
    using namespace OpenSubdiv;

    // convert to mesh
    FarPatchTables::PatchArrayVector const &patchArrays =
        patchTables->GetPatchArrayVector();
    FarPatchTables::PatchParamTable const &patchParam =
        patchTables->GetPatchParamTable();

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
    _mesh.bezierVertices.clear();
    _mesh.bezierBounds.clear();
    _mesh.colors.clear();

    std::vector<int> cpIndices; // 16 * numPatches
    // iterate patch types.
    for (FarPatchTables::PatchArrayVector::const_iterator it = patchArrays.begin();
         it != patchArrays.end(); ++it) {

        int numPatches = 0;
        FarPatchTables::Descriptor desc = it->GetDescriptor();
        VERBOSE("TransitionType = %d\n", desc.GetPattern());

        switch(desc.GetType()) {
        case FarPatchTables::REGULAR:
            numPatches = convertRegular(_mesh.bezierVertices,
                                        _mesh.bezierBounds,
                                        cpIndices,
                                        &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::BOUNDARY:
            numPatches = convertBoundary(_mesh.bezierVertices,
                                         _mesh.bezierBounds,
                                         cpIndices,
                                         &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::CORNER:
            numPatches = convertCorner(_mesh.bezierVertices,
                                       _mesh.bezierBounds,
                                       cpIndices,
                                       &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::GREGORY:
            numPatches = convertGregory(_mesh.bezierVertices,
                                        _mesh.bezierBounds,
                                        cpIndices,
                                        &vertices[0], patchTables, *it);
            break;
        case FarPatchTables::GREGORY_BOUNDARY:
            numPatches = convertBoundaryGregory(_mesh.bezierVertices,
                                                _mesh.bezierBounds,
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
                _mesh.colors.push_back(color[0]);
                _mesh.colors.push_back(color[1]);
                _mesh.colors.push_back(color[2]);
            }
        }
    }

    // for (int i = 0; i < hbrMesh->GetNumFaces(); ++i) {
    //     OsdHbrFace *f = hbrMesh->GetFace(i);
    //     bool out = ((f->isTransitionPatch() or f->_adaptiveFlags.patchType!=OsdHbrFace::kEnd) and
    //                 (not f->_adaptiveFlags.isExtraordinary) and
    //                 (f->_adaptiveFlags.bverts!=1));
    //     VERBOSE("Patch %d out=%x\n", i, out);
    // }


    // vertex position verification pass
    if (_watertight) {
        using namespace OsdUtil;

        std::vector<int> filled(vertices.size());
        std::vector<vec3f> sharedPositions(vertices.size());
        std::vector<int> levels(vertices.size());
        std::map<Edge, Bezier > edgeBeziers;

        for (int i = 0; i < numTotalPatches; ++i) {
            int hbrFace = patchParam[i].hbrFaceIndex;
            OsdHbrFace *face = hbrMesh->GetFace(hbrFace);
            VERBOSE("\n============ patch %d (HBR : %d)==============\n", i, face->GetID());
            int cornerQuad[] = { 0, 3, 12, 15 };
            int parentQuad[] = { 0, 2, 1, 3 };

            int edgeVerts[][4] = { {0, 1, 2, 3}, {0, 4, 8, 12}, {3, 7, 11, 15}, {12, 13, 14, 15} };
//            int edgeVerts[][4] = { {0, 4, 8, 12},{12, 13, 14, 15}, {15, 11, 7, 3}, {3, 2, 1, 0} };

            int edgeParents[][2] = { {0, 2}, {0, 1}, {2, 3}, {1, 3} };
            for (int j = 0; j < 4; ++j) {
                float *v = &_mesh.bezierVertices[(i*16 + cornerQuad[j])*3];
                int parent = cpIndices[i*4+parentQuad[j]];
                if (parent < 0) continue;
                // parent mapping
                // parent = vertexParentIDs[parent];
                parent = farToHbr[parent];
                //VERBOSE("%d/%d  %d: %f %f %f\n", i, cornerQuad[j], parent, v[0], v[1], v[2]);
#if 0
                VERBOSE("%d\n", parent);
                if (filled[parent] == 0) {
                    filled[parent] = 1;
                    sharedPositions[parent] = vec3f(v[0], v[1], v[2]);;
                    levels[parent] = patchParam[i].bitField.GetDepth();
                } else {
                    vec3f pos(v[0], v[1], v[2]);
                    if (pos != sharedPositions[parent]) {
                        VERBOSE("%d/%d ParentID = %d: (level=%d/%d) delta = %g %g %g, (%g, %g, %g)\n",
                               i, cornerQuad[j], parent,
                               patchParam[i].bitField.GetDepth(), levels[parent],
                               sharedPositions[parent][0]-pos[0],
                               sharedPositions[parent][1]-pos[1],
                               sharedPositions[parent][2]-pos[2],
                               pos[0], pos[1], pos[2]);
                        // update points
                        v[0] = sharedPositions[parent][0];
                        v[1] = sharedPositions[parent][1];
                        v[2] = sharedPositions[parent][2];
                    }
                }
#endif

                // for edge verts
                vec3f e0(&_mesh.bezierVertices[(i*16 + edgeVerts[j][0]) * 3]);
                vec3f e1(&_mesh.bezierVertices[(i*16 + edgeVerts[j][1]) * 3]);
                vec3f e2(&_mesh.bezierVertices[(i*16 + edgeVerts[j][2]) * 3]);
                vec3f e3(&_mesh.bezierVertices[(i*16 + edgeVerts[j][3]) * 3]);

                int ep0 = cpIndices[i*4 + edgeParents[j][0]];
                int ep1 = cpIndices[i*4 + edgeParents[j][1]];
                //ep0 = vertexParentIDs[ep0];
                //ep1 = vertexParentIDs[ep1];
                ep0 = farToHbr[ep0];
                ep1 = farToHbr[ep1];
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
                } else {
#if 0
                    std::pair<vec3f, vec3f> ref = edgePositions[edge];
                    vec3f re0 = ref.first;
                    vec3f re1 = ref.second;
                    if (edge._edges[0] != ep0) std::swap(re0, re1);

                    if (e0 != re0) {
                        VERBOSE("Edge %d-%d: (level=%d) delta = %g %g %g\n",
                               ep0, ep1,
                               patchParam[i].bitField.GetDepth(),
                               re0[0] - e0[0],
                               re0[1] - e0[1],
                               re0[2] - e0[2]);
                        ppe0[0] = re0[0];
                        ppe0[1] = re0[1];
                        ppe0[2] = re0[2];
                    }
                    if (e1 != re1) {
                        VERBOSE("Edge %d-%d: (level=%d) delta = %g %g %g\n",
                               ep0, ep1,
                               patchParam[i].bitField.GetDepth(),
                               re1[0] - e1[0],
                               re1[1] - e1[1],
                               re1[2] - e1[2]);
                        ppe1[0] = re1[0];
                        ppe1[1] = re1[1];
                        ppe1[2] = re1[2];
                    }
#endif
                }
            }
        }

        // different level subdivision
        // water-tight critical
        for (int i = 0; i < numTotalPatches; ++i) {
            int hbrFace = patchParam[i].hbrFaceIndex;
            OsdHbrFace *face = hbrMesh->GetFace(hbrFace);

            if (not face->_adaptiveFlags.isCritical) continue;

            VERBOSE("Critical patch %d, %d----\n", i, hbrFace);

            OsdHbrFace *parentFace = face->GetParent();

            if (not parentFace) continue;

            if (parentFace->GetNumVertices() != 4) {
                VERBOSE("non-quad parent %d\n", parentFace->GetID());
                continue;
            }
            // printf("parent face = <%d, lv=%d>, extraordinary = %d ", parentFace->GetID(), parentFace->GetDepth(), parentFace->_adaptiveFlags.isExtraordinary);
            // for (int j = 0; j < 4; ++j) {
            //     printf("%d, ", parentFace->GetVertex(j)->GetID());
            // }
            // printf("\n");

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

            // for (int j = 0; j < 4; ++j) {
            //     int edgeParents[][2] = { {0, 2}, {0, 1}, {2, 3}, {1, 3} };
            //     int ep0 = cpIndices[i*4 + edgeParents[j][0]];
            //     int ep1 = cpIndices[i*4 + edgeParents[j][1]];
            //     ep0 = farToHbr[ep0];
            //     ep1 = farToHbr[ep1];
            //     VERBOSE("Edge %d-%d\n", ep0, ep1);
            // }

            // childindex 0 -> edge 3, 0
            // childindex 1 -> edge 0, 1
            // childindex 2 -> edge 1, 2
            // childindex 3 -> edge 2, 3
            // for each edge
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

                if (edgeBeziers.count(e) != 0) {
                    Bezier parentEdge = edgeBeziers[e];  // todo, lookup
                    if (e._edges[0] != v0->GetID()) {
                        parentEdge.Reverse();
                    }

                    // find corresponding edge at patch[i]
                    // cut original bezier
                    vec3f tmp[2][4];
                    OsdUtil::BezierSplit<4, 0, /*stride*/1, vec3f, float> split(
                        tmp[0], tmp[1], edgeBeziers[e].cp, 0.5);

                    VERBOSE("[[%f %f %f - ", tmp[0][0][0], tmp[0][0][1], tmp[0][0][2]);
                    VERBOSE("%f %f %f - ", tmp[1][0][0], tmp[1][0][1], tmp[1][0][2]);
                    VERBOSE("%f %f %f]]\n", tmp[1][3][0], tmp[1][3][1], tmp[1][3][2]);

                    if (consolidateBezier(Bezier(tmp[0]), Bezier(tmp[1]),
                                          &_mesh.bezierVertices[(i*16*3)]) == false) {
                        printf("Matching error. face = %d, hbr face=%d\n", i, face->GetID());
                    }
                } else {
                    printf("Topology error --- Not found hbr verts in edge dict%d,%d\n",
                           v0->GetID(), v1->GetID());
                }
            }
        }
    }


    _mesh.numTriangles = 0;
    _mesh.numBezierPatches = numTotalPatches;
    _mesh.patchParams = &(patchParam[0]);

    assert(numTotalPatches*16*3 == (int)_mesh.bezierVertices.size());
}

static void evalBezier(float *p, float *n, float u, float v, const float *cp)
{
    OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3f, float, 4> patch((const OsdUtil::vec3f*)cp);
    OsdUtil::vec3f b = patch.Evaluate(u, v);
    p[0] = b[0];
    p[1] = b[1];
    p[2] = b[2];
    OsdUtil::vec3f du = patch.EvaluateDu(u, v);
    OsdUtil::vec3f dv = patch.EvaluateDv(u, v);
    n[0] = -(du[1] * dv[2] - du[2] * dv[1]);
    n[1] = -(du[2] * dv[0] - du[0] * dv[2]);
    n[2] = -(du[0] * dv[1] - du[1] * dv[0]);
}

void
Scene::Tessellate(int level)
{
    float *cp = &_mesh.bezierVertices[0];
    int vindex = 0;
    int numTriangles = 0;

    _mesh.triVertices.clear();
    _mesh.triNormals.clear();
    _mesh.faces.clear();
    _mesh.colors.clear();
    std::vector<float> colors;
    for (size_t patchIndex = 0; patchIndex < _mesh.numBezierPatches; ++patchIndex) {

        const OpenSubdiv::FarPatchParam &param = _mesh.patchParams[patchIndex];
        unsigned int bits = param.bitField.field;
        int patchLevel = (bits & 0xf);
        int div = 1 << std::max(0, level-patchLevel);

        float p[3], n[3];
        for (int u = 0; u <= div; ++u) {
            for (int v = 0; v <= div; ++v) {
                evalBezier(p, n, u/float(div), v/float(div), cp);
                _mesh.triVertices.push_back(p[0]);
                _mesh.triVertices.push_back(p[1]);
                _mesh.triVertices.push_back(p[2]);
                _mesh.triNormals.push_back(n[0]);
                _mesh.triNormals.push_back(n[1]);
                _mesh.triNormals.push_back(n[2]);

                if (u < div and v < div) {
                    _mesh.faces.push_back(vindex);
                    _mesh.faces.push_back(vindex+1);
                    _mesh.faces.push_back(vindex+div+1);

                    _mesh.faces.push_back(vindex+div+1);
                    _mesh.faces.push_back(vindex+1);
                    _mesh.faces.push_back(vindex+div+2);

                    colors.push_back(_mesh.colors[patchIndex*3+0]);
                    colors.push_back(_mesh.colors[patchIndex*3+1]);
                    colors.push_back(_mesh.colors[patchIndex*3+2]);
                    colors.push_back(_mesh.colors[patchIndex*3+0]);
                    colors.push_back(_mesh.colors[patchIndex*3+1]);
                    colors.push_back(_mesh.colors[patchIndex*3+2]);
                    numTriangles += 2;
                }
                ++vindex;
            }
        }

        // next
        cp += 16*3;
    }

    _mesh.numTriangles = numTriangles;
    _mesh.colors.swap(colors);

    // clear bezier patch
    _mesh.numBezierPatches = 0;
    _mesh.bezierVertices.clear();
    _mesh.bezierBounds.clear();
}

void
Scene::Build()
{
    BVHBuildOptions options; // Use default option

    printf("  BVH build option:\n");
    printf("    # of leaf primitives: %d\n", options.minLeafPrimitives);
    printf("    SAH binsize         : %d\n", options.binSize);

    printf("  # of triangles : %ld\n", _mesh.numTriangles);
    printf("  # of bezier patches : %ld\n", _mesh.numBezierPatches);

    _accel = BVHAccel();
    _accel.Build(&_mesh, options);

    BVHBuildStatistics stats = _accel.GetStatistics();

    printf("  BVH statistics:\n");
    printf("    # of leaf   nodes: %d\n", stats.numLeafNodes);
    printf("    # of branch nodes: %d\n", stats.numBranchNodes);
    printf("  Max tree depth   : %d\n", stats.maxTreeDepth);
}

void
Scene::VBOBuild()
{
    // create vbo for display
    if (_vbo == 0) glGenBuffers(1, &_vbo);
    glBindBuffer(GL_ARRAY_BUFFER, _vbo);
    glBufferData(GL_ARRAY_BUFFER, _accel.GetNodes().size() * sizeof(BVHNode),
                 &_accel.GetNodes()[0], GL_STATIC_DRAW);
    glBindBuffer(GL_ARRAY_BUFFER, 0);
}

#ifdef OPENSUBDIV_HAS_TBB
class Kernel {
public:
    Kernel(int width, int stepIndex, int step, BVHAccel *accel, Mesh *mesh,
           Camera *camera, float *image, Scene *scene) :
        _width(width), _stepIndex(stepIndex), _step(step), _accel(accel),
        _mesh(mesh), _camera(camera), _image(image), _scene(scene) {
    }

    void operator() (tbb::blocked_range<int> const &r) const {
        for (int rr = r.begin(); rr < r.end(); ++rr) {
            int y = rr*_step + _stepIndex/_step;
            for (int x = _stepIndex%_step; x < _width; x += _step) {

                float u = 0.5;
                float v = 0.5;

                Ray ray = _camera->GenerateRay(x + u + _step / 2.0f, y + v + _step / 2.0f);

                Intersection isect;
                bool hit = _accel->Traverse(isect, _mesh, ray);

                float rgba[4] = { 0, 0, 0, 0 };
                float *d = _image + 4 * (y * _width + x);
                if (hit) {
                    _scene->Shade(rgba, isect, ray);
                } else {
                    // Maya like gradation. Maybe helpful to check crack visually.
                    rgba[0] = 0.1f; 
                    rgba[1] = 0.1f;
                    rgba[2] = 0.4f * ((_width - x - 1)/(double)_width);
                    rgba[3] = 1.0f;
                }
                d[0] = rgba[0];
                d[1] = rgba[1];
                d[2] = rgba[2];
                d[3] = rgba[3];
            }
        }
    }

private:
    int _width;
    int _stepIndex;
    int _step;
    BVHAccel *_accel;
    Mesh *_mesh;
    Camera *_camera;
    float *_image;
    Scene *_scene;
};
#endif

void
Scene::Render(int width, int height, double fov,
              std::vector<float> &image, // RGB
              const float eye[3],
              const float lookat[3], const float up[3],
              int step, int stepIndex, int intersectKernel, float uvMargin)
{
    _accel.SetIntersectKernel(intersectKernel);
    _accel.SetUVMargin(uvMargin);

    Camera camera;

    double deye[3] = { eye[0], eye[1], eye[2] };
    double dlook[3] = { lookat[0], lookat[1], lookat[2] };
    double dup[3] = { up[0], up[1], up[2] };

    camera.BuildCameraFrame(deye, dlook, dup, fov, width, height);
    assert((int)image.size() >= 3 * width * height);

#ifdef OPENSUBDIV_HAS_TBB
    tbb::blocked_range<int> range(0, height/step, 1);

    Kernel kernel(width, stepIndex, step, &_accel, &_mesh, &camera, &image[0], this);
    tbb::parallel_for(range, kernel);
#else

#ifdef _OPENMP
#pragma omp parallel for schedule(dynamic, 1)
#endif
    for (int y = stepIndex/step; y < height; y += step) {

        for (int x = stepIndex%step; x < width; x += step) {

            float u = 0.5;
            float v = 0.5;

            Ray ray = camera.GenerateRay(x + u + step / 2.0f, y + v + step / 2.0f);

            Intersection isect;
            bool hit = _accel.Traverse(isect, &_mesh, ray);

            float rgba[4] = { 0, 0, 0, 0};
            if (hit) {
                Shade(rgba, isect, ray);
            } else {
                // Maya like gradation. Maybe helpful to check crack visually.
                rgba[0] = 0.1f; 
                rgba[1] = 0.1f;
                rgba[2] = 0.4f * ((height - y - 1)/(double)height);
                rgba[3] = 1.0f;
            }
            image[4 * (y * width + x) + 0] = rgba[0];
            image[4 * (y * width + x) + 1] = rgba[1];
            image[4 * (y * width + x) + 2] = rgba[2];
            image[4 * (y * width + x) + 3] = rgba[3];
        }
    }
#endif
}

inline float randomreal(void) {
  static unsigned int x = 123456789, y = 362436069, z = 521288629,
      w = 88675123;
  unsigned t = x ^ (x << 11);
  x = y;
  y = z;
  z = w;
  w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
  return w * (1.0f / 4294967296.0f);
}

void
Scene::Shade(float rgba[4], const Intersection &isect, const Ray &ray)
{
    real3 I = ray.dir;

    real d = std::max(real(0), vdot(I, isect.normal));
    real3 color;
    if (_mode == SHADED) {
        //        real3 reflect = I - 2 * d * isect.normal;
        real s = 0;//pow(std::max(0.0f, -vdot(ray.dir, reflect)), 32);
        color = d * real3(0.8, 0.8, 0.8) + s * real3(1, 1, 1);
    } else if (_mode == PTEX_COORD) {
        color = real3(isect.u, isect.v, 1);
    } else if (_mode == PATCH_TYPE) {
        float l = isect.level * 0.05;
        color = d * (real3(&_mesh.colors[isect.patchID*3])
                     - real3(l, l, l));
        color[0] = std::max(real(0), color[0]);
        color[1] = std::max(real(0), color[1]);
        color[2] = std::max(real(0), color[2]);
    } else if (_mode == CLIP_LEVEL) {
        float l = isect.clipLevel * 0.05;
        color = (real3(&_mesh.colors[isect.patchID*3])
                 - real3(l, l, l));
        color[0] = std::max(real(0), color[0]);
        color[1] = std::max(real(0), color[1]);
        color[2] = std::max(real(0), color[2]);
    } else if (_mode == AO) {
        Intersection si;
        Ray sray;
        int numHits = 0;

        int numSamples = 16;
        for (int i = 0; i < numSamples; ++i) {
            real3 sample = real3(0.5-randomreal(), 0.5-randomreal(), 0.5-randomreal());
            sample.normalize();
            sray.dir = sample * (vdot(sample, isect.normal) > 0 ? -1 : 1);
            sray.invDir = sray.dir.neg();
            sray.org = ray.org + ray.dir * isect.t + sray.dir * 0.0001;
            numHits += _accel.Traverse(si, &_mesh, sray) ? 1 : 0;
        }

        color[0] = color[1] = color[2] = d * (1.0-numHits/float(numSamples));
    } else if (_mode == TRANSPARENT) {
        float alpha = 0.25 * (1.0 - rgba[3]);
        rgba[0] += d * alpha;
        rgba[1] += d * alpha;
        rgba[2] += d * alpha;
        rgba[3] += alpha;

        if (alpha < 0.9) {
            Intersection si;
            Ray sray;
            sray.dir = ray.dir;
            sray.invDir = ray.invDir;
            sray.org = ray.org + ray.dir * (isect.t + 0.0001);
            if (_accel.Traverse(si, &_mesh, sray)) {
                Shade(rgba, si, sray);
            } else {
                rgba[3] = 1.0;
                return;
            }
        } else {
            rgba[3] = 1.0;
        }
        return;
    }
    rgba[0] = color[0];
    rgba[1] = color[1];
    rgba[2] = color[2];
    rgba[3] = 1.0;
}

