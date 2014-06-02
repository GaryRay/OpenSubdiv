//
//   Copyright 2013 Pixar
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

#include <stdio.h>
#include <float.h>

#include <far/meshFactory.h>
#include <far/dispatcher.h>

#include "../regression/common/shape_utils.h"
#include "eson.h"
#include "bezier.h"

//
// Regression testing matching Far to Hbr (default CPU implementation)
//
// Notes:
// - precision is currently held at 1e-6
//
// - results cannot be bitwise identical as some vertex interpolations
//   are not happening in the same order.
//
// - only vertex interpolation is being tested at the moment.
//
#define PRECISION 1e-6

//------------------------------------------------------------------------------
// Vertex class implementation
struct xyzVV {

    xyzVV() { }

    xyzVV( int /*i*/ ) { }

    xyzVV( float x, float z, float y ) { _pos[0]=x; _pos[1]=y; _pos[2]=z; }

    xyzVV( const xyzVV & src ) { _pos[0]=src._pos[0]; _pos[1]=src._pos[1]; _pos[2]=src._pos[2]; }

   ~xyzVV( ) { }

    void AddWithWeight(const xyzVV& src, float weight) { 
        _pos[0]+=weight*src._pos[0]; 
        _pos[1]+=weight*src._pos[1]; 
        _pos[2]+=weight*src._pos[2]; 
    }

    void AddVaryingWithWeight(const xyzVV& , float) { }

    void Clear( void * =0 ) { _pos[0]=_pos[1]=_pos[2]=0.0f; }

    void SetPosition(float x, float z, float y) { _pos[0]=x; _pos[1]=y; _pos[2]=z; }

    void ApplyVertexEdit(const OpenSubdiv::HbrVertexEdit<xyzVV> & edit) {
        const float *src = edit.GetEdit();
        switch(edit.GetOperation()) {
          case OpenSubdiv::HbrHierarchicalEdit<xyzVV>::Set:
            _pos[0] = src[0];
            _pos[1] = src[1];
            _pos[2] = src[2];
            break;
          case OpenSubdiv::HbrHierarchicalEdit<xyzVV>::Add:
            _pos[0] += src[0];
            _pos[1] += src[1];
            _pos[2] += src[2];
            break;
          case OpenSubdiv::HbrHierarchicalEdit<xyzVV>::Subtract:
            _pos[0] -= src[0];
            _pos[1] -= src[1];
            _pos[2] -= src[2];
            break;
        }
    }

    void ApplyVertexEdit(OpenSubdiv::FarVertexEdit const & edit) {
        const float *src = edit.GetEdit();
        switch(edit.GetOperation()) {
          case OpenSubdiv::FarVertexEdit::Set:
            _pos[0] = src[0];
            _pos[1] = src[1];
            _pos[2] = src[2];
            break;
          case OpenSubdiv::FarVertexEdit::Add:
            _pos[0] += src[0];
            _pos[1] += src[1];
            _pos[2] += src[2];
            break;
        }
    }
    
    void ApplyMovingVertexEdit(const OpenSubdiv::HbrMovingVertexEdit<xyzVV> &) { }

    const float * GetPos() const { return _pos; }

private:
    float _pos[3];
};

//------------------------------------------------------------------------------
class xyzFV;
typedef OpenSubdiv::HbrMesh<xyzVV>           xyzmesh;
typedef OpenSubdiv::HbrFace<xyzVV>           xyzface;
typedef OpenSubdiv::HbrVertex<xyzVV>         xyzvertex;
typedef OpenSubdiv::HbrHalfedge<xyzVV>       xyzhalfedge;
typedef OpenSubdiv::HbrFaceOperator<xyzVV>   xyzFaceOperator;
typedef OpenSubdiv::HbrVertexOperator<xyzVV> xyzVertexOperator;

typedef OpenSubdiv::FarMesh<xyzVV>              fMesh;
typedef OpenSubdiv::FarMeshFactory<xyzVV>       fMeshFactory;
typedef OpenSubdiv::FarSubdivisionTables        fSubdivision;
typedef OpenSubdiv::FarPatchTables              fPatches;

static int g_level = 4;


//------------------------------------------------------------------------------
int checkMesh( char const * name, xyzmesh * hmesh, int levels, Scheme scheme=kCatmark ) {

    using namespace OpenSubdiv;

    assert(name);

    int count=0;

    fMeshFactory fact( hmesh, levels, /*adaptive=*/true);
    fMesh * m = fact.Create( );
    static OpenSubdiv::FarComputeController computeController;
    computeController.Refine(m);

    // dump patch as eson

    // centering vertices.
    std::vector<float> vertices;
    {
        float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
        float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};
        for (int i = 0; i < (int)m->GetVertices().size(); ++i) {
            xyzVV const &v = m->GetVertices()[i];
            for (int j = 0; j < 3; ++j) {
                min[j] = std::min(min[j], v.GetPos()[j]);
                max[j] = std::max(max[j], v.GetPos()[j]);
            }
        }
        float center[3] = { (max[0]+min[0])*0.5,
                            (max[1]+min[1])*0.5,
                            (max[2]+min[2])*0.5 };
        for (int i = 0; i < (int)m->GetVertices().size(); ++i) {
            xyzVV v = m->GetVertices()[i];
            float x = v.GetPos()[0];
            float y = v.GetPos()[1];
            float z = v.GetPos()[2];
            vertices.push_back(x-center[0]);
            vertices.push_back(y-center[1]);
            vertices.push_back(z-center[2]);
        }
    }

    FarPatchTables const *patchTables = m->GetPatchTables();
    FarPatchTables::PatchArrayVector const &patchArrays = patchTables->GetPatchArrayVector();
    FarPatchTables::PatchParamTable const &patchParam = patchTables->GetPatchParamTable();

    int numPatches = 0;
    std::vector<float> bezierVertices;
    // iterate patch types.
    printf("begin\n");
    for (FarPatchTables::PatchArrayVector::const_iterator it = patchArrays.begin();
         it != patchArrays.end(); ++it) {

        switch(it->GetDescriptor().GetType()) {
        case FarPatchTables::REGULAR:
            numPatches += convertRegular(bezierVertices, vertices, patchTables, *it);
            break;
        case FarPatchTables::BOUNDARY:
            numPatches += convertBoundary(bezierVertices, vertices, patchTables, *it);
            break;
        case FarPatchTables::CORNER:
            numPatches += convertCorner(bezierVertices, vertices, patchTables, *it);
            break;
        case FarPatchTables::GREGORY:
            numPatches += convertGregory(bezierVertices, vertices, patchTables, *it);
            break;
        case FarPatchTables::GREGORY_BOUNDARY:
            numPatches += convertBoundaryGregory(bezierVertices, vertices, patchTables, *it);
            break;
        default:
            break;
        }

    }

    eson::Object mesh;
    int nverts = m->GetNumVertices();
    mesh["num_vertices"] =
        eson::Value((int64_t)nverts);
    mesh["vertices"] =
        eson::Value((uint8_t*)&vertices[0], sizeof(float)*nverts*3);
    mesh["num_bezier_patches"] =
        eson::Value((int64_t)numPatches);
    mesh["patch_param"] =
        eson::Value((uint8_t*)&patchParam[0], sizeof(OpenSubdiv::FarPatchParam)*patchParam.size());
    mesh["bezier_vertices"] =
        eson::Value((uint8_t*)&bezierVertices[0], sizeof(float)*bezierVertices.size());

    printf("%s, verts=%d, patches=%d, %d\n", name, nverts, numPatches, (int)bezierVertices.size());

    assert(numPatches*16*3 == (int)bezierVertices.size());

    eson::Value v = eson::Value(mesh);
    int64_t size = v.Size();

    std::vector<uint8_t> buf(size);
    uint8_t* ptr = &buf[0]; 

    ptr = v.Serialize(ptr);
    assert((ptr-&buf[0]) == size);

    std::string filename = std::string(name) + ".eson";
    FILE* fp = fopen(filename.c_str(), "wb");
    fwrite(&buf[0], 1, size, fp);
    fclose(fp);

    delete hmesh;
    delete m;

    return count;
}

//------------------------------------------------------------------------------
static void parseArgs(int argc, char ** argv) {
    if (argc>1) {
        for (int i=1; i<argc; ++i) {
            if (strcmp(argv[i],"-l")==0) {
                g_level = atoi(argv[++i]);
            } else {
                printf("Unknown argument \"%s\".\n", argv[i]);
                exit(1);
            }
        }
    }
}

//------------------------------------------------------------------------------
int main(int argc, char ** argv) {

    int total=0;

    parseArgs(argc, argv);

    int levels = g_level;

#define test_catmark_edgeonly
#define test_catmark_edgecorner
#define test_catmark_flap
#define test_catmark_pyramid
#define test_catmark_pyramid_creases0
#define test_catmark_pyramid_creases1
#define test_catmark_cube
#define test_catmark_cube_creases0
#define test_catmark_cube_creases1
#define test_catmark_cube_corner0
#define test_catmark_cube_corner1
#define test_catmark_cube_corner2
#define test_catmark_cube_corner3
#define test_catmark_cube_corner4
#define test_catmark_dart_edgeonly
#define test_catmark_dart_edgecorner
#define test_catmark_tent
#define test_catmark_tent_creases0
#define test_catmark_tent_creases1
#define test_catmark_square_hedit0
#define test_catmark_square_hedit1
#define test_catmark_square_hedit2
#define test_catmark_square_hedit3
#define test_catmark_car
#define test_catmark_bishop
#define test_catmark_pawn
#define test_catmark_rook
#define test_catmark_helmet
#define test_catmark_gregory_test1
#define test_catmark_gregory_test2
#define test_catmark_gregory_test3
#define test_catmark_gregory_test4

#if 0
#define test_loop_triangle_edgeonly
#define test_loop_triangle_edgecorner
#define test_loop_icosahedron
#define test_loop_cube
#define test_loop_cube_creases0
#define test_loop_cube_creases1

#define test_bilinear_cube
#endif

#ifdef test_catmark_edgeonly
#include "../regression/shapes/catmark_edgeonly.h"
    total += checkMesh( "test_catmark_edgeonly", simpleHbr<xyzVV>(catmark_edgeonly.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_edgecorner
#include "../regression/shapes/catmark_edgecorner.h"
    total += checkMesh( "test_catmark_edgecorner", simpleHbr<xyzVV>(catmark_edgecorner.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_pyramid
#include "../regression/shapes/catmark_pyramid.h"
    total += checkMesh( "test_catmark_pyramid", simpleHbr<xyzVV>(catmark_pyramid.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_pyramid_creases0
#include "../regression/shapes/catmark_pyramid_creases0.h"
    total += checkMesh( "test_catmark_pyramid_creases0", simpleHbr<xyzVV>(catmark_pyramid_creases0.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_pyramid_creases1
#include "../regression/shapes/catmark_pyramid_creases1.h"
    total += checkMesh( "test_catmark_pyramid_creases1", simpleHbr<xyzVV>(catmark_pyramid_creases1.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube
#include "../regression/shapes/catmark_cube.h"
    total += checkMesh( "test_catmark_cube", simpleHbr<xyzVV>(catmark_cube.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_creases0
#include "../regression/shapes/catmark_cube_creases0.h"
    total += checkMesh( "test_catmark_cube_creases0", simpleHbr<xyzVV>(catmark_cube_creases0.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_creases1
#include "../regression/shapes/catmark_cube_creases1.h"
    total += checkMesh( "test_catmark_cube_creases1", simpleHbr<xyzVV>(catmark_cube_creases1.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_corner0
#include "../regression/shapes/catmark_cube_corner0.h"
    total += checkMesh( "test_catmark_cube_corner0", simpleHbr<xyzVV>(catmark_cube_corner0.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_corner1
#include "../regression/shapes/catmark_cube_corner1.h"
    total += checkMesh( "test_catmark_cube_corner1", simpleHbr<xyzVV>(catmark_cube_corner1.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_corner2
#include "../regression/shapes/catmark_cube_corner2.h"
    total += checkMesh( "test_catmark_cube_corner2", simpleHbr<xyzVV>(catmark_cube_corner2.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_corner3
#include "../regression/shapes/catmark_cube_corner3.h"
    total += checkMesh( "test_catmark_cube_corner3", simpleHbr<xyzVV>(catmark_cube_corner3.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_cube_corner4
#include "../regression/shapes/catmark_cube_corner4.h"
    total += checkMesh( "test_catmark_cube_corner4", simpleHbr<xyzVV>(catmark_cube_corner4.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_dart_edgecorner
#include "../regression/shapes/catmark_dart_edgecorner.h"
    total += checkMesh( "test_catmark_dart_edgecorner", simpleHbr<xyzVV>(catmark_dart_edgecorner.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_dart_edgeonly
#include "../regression/shapes/catmark_dart_edgeonly.h"
    total += checkMesh( "test_catmark_dart_edgeonly", simpleHbr<xyzVV>(catmark_dart_edgeonly.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_flap
#include "../regression/shapes/catmark_flap.h"
    total += checkMesh( "test_catmark_flap", simpleHbr<xyzVV>(catmark_flap.c_str(), kCatmark, 0), levels);
#endif

#ifdef test_catmark_tent
#include "../regression/shapes/catmark_tent.h"
    total += checkMesh( "test_catmark_tent", simpleHbr<xyzVV>(catmark_tent.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_tent_creases0
#include "../regression/shapes/catmark_tent_creases0.h"
    total += checkMesh( "test_catmark_tent_creases0", simpleHbr<xyzVV>(catmark_tent_creases0.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_tent_creases1
#include "../regression/shapes/catmark_tent_creases1.h"
    total += checkMesh( "test_catmark_tent_creases1", simpleHbr<xyzVV>(catmark_tent_creases1.c_str(), kCatmark, NULL), levels );
#endif

#ifdef test_catmark_square_hedit0
#include "../regression/shapes/catmark_square_hedit0.h"
    total += checkMesh( "test_catmark_square_hedit0", simpleHbr<xyzVV>(catmark_square_hedit0.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_square_hedit1
#include "../regression/shapes/catmark_square_hedit1.h"
    total += checkMesh( "test_catmark_square_hedit1", simpleHbr<xyzVV>(catmark_square_hedit1.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_square_hedit2
#include "../regression/shapes/catmark_square_hedit2.h"
    total += checkMesh( "test_catmark_square_hedit2", simpleHbr<xyzVV>(catmark_square_hedit2.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_square_hedit3
#include "../regression/shapes/catmark_square_hedit3.h"
    total += checkMesh( "test_catmark_square_hedit3", simpleHbr<xyzVV>(catmark_square_hedit3.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_car
#include "../regression/shapes/catmark_car.h"
    total += checkMesh( "test_catmark_car", simpleHbr<xyzVV>(catmark_car.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_bishop
#include "../regression/shapes/catmark_bishop.h"
    total += checkMesh( "test_catmark_bishop", simpleHbr<xyzVV>(catmark_bishop.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_helmet
#include "../regression/shapes/catmark_helmet.h"
    total += checkMesh( "test_catmark_helmet", simpleHbr<xyzVV>(catmark_helmet.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_pawn
#include "../regression/shapes/catmark_pawn.h"
    total += checkMesh( "test_catmark_pawn", simpleHbr<xyzVV>(catmark_pawn.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_rook
#include "../regression/shapes/catmark_rook.h"
    total += checkMesh( "test_catmark_rook", simpleHbr<xyzVV>(catmark_rook.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_gregory_test1
#include "../regression/shapes/catmark_gregory_test1.h"
    total += checkMesh( "test_catmark_gregory_test1", simpleHbr<xyzVV>(catmark_gregory_test1.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_gregory_test2
#include "../regression/shapes/catmark_gregory_test2.h"
    total += checkMesh( "test_catmark_gregory_test2", simpleHbr<xyzVV>(catmark_gregory_test2.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_gregory_test3
#include "../regression/shapes/catmark_gregory_test3.h"
    total += checkMesh( "test_catmark_gregory_test3", simpleHbr<xyzVV>(catmark_gregory_test3.c_str(), kCatmark, 0), levels );
#endif

#ifdef test_catmark_gregory_test4
#include "../regression/shapes/catmark_gregory_test4.h"
    total += checkMesh( "test_catmark_gregory_test4", simpleHbr<xyzVV>(catmark_gregory_test4.c_str(), kCatmark, 0), levels );
#endif


#ifdef test_loop_triangle_edgeonly
#include "../regression/shapes/loop_triangle_edgeonly.h"
    total += checkMesh( "test_loop_triangle_edgeonly", simpleHbr<xyzVV>(loop_triangle_edgeonly.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_triangle_edgecorner
#include "../regression/shapes/loop_triangle_edgecorner.h"
    total += checkMesh( "test_loop_triangle_edgecorner", simpleHbr<xyzVV>(loop_triangle_edgecorner.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_saddle_edgeonly
#include "../regression/shapes/loop_saddle_edgeonly.h"
    total += checkMesh( "test_loop_saddle_edgeonly", simpleHbr<xyzVV>(loop_saddle_edgeonly.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_saddle_edgecorner
#include "../regression/shapes/loop_saddle_edgecorner.h"
    total += checkMesh( "test_loop_saddle_edgecorner", simpleHbr<xyzVV>(loop_saddle_edgecorner.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_icosahedron
#include "../regression/shapes/loop_icosahedron.h"
    total += checkMesh( "test_loop_icosahedron", simpleHbr<xyzVV>(loop_icosahedron.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_cube
#include "../regression/shapes/loop_cube.h"
    total += checkMesh( "test_loop_cube", simpleHbr<xyzVV>(loop_cube.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_cube_creases0
#include "../regression/shapes/loop_cube_creases0.h"
    total += checkMesh( "test_loop_cube_creases0", simpleHbr<xyzVV>(loop_cube_creases0.c_str(), kLoop, 0), levels, kLoop );
#endif

#ifdef test_loop_cube_creases1
#include "../regression/shapes/loop_cube_creases1.h"
    total += checkMesh( "test_loop_cube_creases1", simpleHbr<xyzVV>(loop_cube_creases1.c_str(), kLoop, 0), levels, kLoop );
#endif


#ifdef test_bilinear_cube
#include "../regression/shapes/bilinear_cube.h"
    total += checkMesh( "test_bilinear_cube", simpleHbr<xyzVV>(bilinear_cube.c_str(), kBilinear, 0), levels, kBilinear );
#endif
}

//------------------------------------------------------------------------------
