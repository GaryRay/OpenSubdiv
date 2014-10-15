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

#include "../osd/tbbSmoothNormalController.h"

#include <math.h>
#include <string.h>
#include <tbb/parallel_for.h>

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Osd {

inline void
cross(float *n, const float *p0, const float *p1, const float *p2) {

    float a[3] = { p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2] };
    float b[3] = { p2[0]-p0[0], p2[1]-p0[1], p2[2]-p0[2] };
    n[0] = a[1]*b[2]-a[2]*b[1];
    n[1] = a[2]*b[0]-a[0]*b[2];
    n[2] = a[0]*b[1]-a[1]*b[0];

    float rn = 1.0f/sqrtf(n[0]*n[0] + n[1]*n[1] + n[2]*n[2]);
    n[0] *= rn;
    n[1] *= rn;
    n[2] *= rn;
}

#define grain_size  200

// TBB kernel to reset normals to 0.0f
class TBBResetKernel {

    float * _oBuffer;
    int _oStride;

public:

    void operator() (tbb::blocked_range<int> const &r) const {

        float * dst = _oBuffer + r.begin() * _oStride;

        // reset normals to 0
        for (int i=r.begin(); i<r.end(); ++i, dst+=_oStride) {
            memset(dst, 0, 3*sizeof(float));
        }
    }

    TBBResetKernel(TBBResetKernel const & other) {
        this->_oBuffer = other._oBuffer;
        this->_oStride = other._oStride;
    }

    TBBResetKernel(float * oBuffer, int oStride) :
        _oBuffer(oBuffer), _oStride(oStride) {
    }
};

// TBB kernel that averages face normals into vertex buffer
class TBBSmoothNormalKernel {

    float const * _iBuffer;
    float * _oBuffer;

    Far::Index const * _vertIndices;

    int _iStride,
        _oStride,
        _numVertices;

public:

    void operator() (tbb::blocked_range<int> const &r) const {

            int idx = r.begin()*_numVertices;

            for (int i=r.begin(); i<r.end(); ++i, idx+=_numVertices) {

                float const * p0 = _iBuffer + _vertIndices[idx+0]*_iStride,
                            * p1 = _iBuffer + _vertIndices[idx+1]*_iStride,
                            * p2 = _iBuffer + _vertIndices[idx+2]*_iStride;

                // compute face normal
                float n[3];
                cross( n, p0, p1, p2 );

                // add normal to all vertices of the face
                for (int j=0; j<_numVertices; ++j) {

                    float * dst = _oBuffer + _vertIndices[idx+j]*_oStride;

                    dst[0] += n[0];
                    dst[1] += n[1];
                    dst[2] += n[2];
                }
            }
    }

    TBBSmoothNormalKernel( TBBSmoothNormalKernel const & other ) {
        this->_iBuffer     = other._iBuffer;
        this->_oBuffer     = other._oBuffer;
        this->_vertIndices = other._vertIndices;
        this->_iStride     = other._iStride;
        this->_oStride     = other._oStride;
        this->_numVertices = other._numVertices;
    }

    TBBSmoothNormalKernel( float const * iBuffer,
                           int iStride,
                           float * oBuffer,
                           int oStride,
                           Far::Index const * vertIndices,
                           int numVertices ) :
        _iBuffer(iBuffer),
        _oBuffer(oBuffer),
        _vertIndices(vertIndices),
        _iStride(iStride),
        _oStride(oStride),
        _numVertices(numVertices) {
    }
};

void TbbSmoothNormalController::_smootheNormals(
    CpuSmoothNormalContext * context) {

    VertexBufferDescriptor const & iDesc = context->GetInputVertexDescriptor(),
                                    & oDesc = context->GetOutputVertexDescriptor();

    assert(iDesc.length==3 and oDesc.length==3);

    float const * iBuffer = context->GetCurrentInputVertexBuffer() + iDesc.offset;
    float * oBuffer = context->GetCurrentOutputVertexBuffer() + oDesc.offset;

    std::vector<Far::Index> const & verts = context->GetControlVertices();

    Far::PatchTables::PatchArrayVector const & parrays = context->GetPatchArrayVector();

    if (verts.empty() or parrays.empty() or (not iBuffer) or (not oBuffer)) {
        return;
    }

    for (int i=0; i<(int)parrays.size(); ++i) {

        Far::PatchTables::PatchArray const & pa = parrays[i];

        Far::PatchTables::Type type = pa.GetDescriptor().GetType();

        if (type==Far::PatchTables::QUADS or type==Far::PatchTables::TRIANGLES) {


            // if necessary, reset all normal values to 0
            if (context->GetResetMemory()) {

                TBBResetKernel resetKernel(oBuffer, oDesc.stride);
                tbb::blocked_range<int> range(0, context->GetNumVertices(), grain_size);
                tbb::parallel_for(range, resetKernel);
            }

            {
                int nv = Far::PatchTables::Descriptor::GetNumControlVertices(type);
                TBBSmoothNormalKernel smoothNormalkernel( iBuffer,
                                                          iDesc.stride,
                                                          oBuffer,
                                                          oDesc.stride,
                                                          &verts[pa.GetVertIndex()],
                                                          nv);

                tbb::blocked_range<int> range(0, pa.GetNumPatches(), grain_size);
                tbb::parallel_for(range, smoothNormalkernel);
            }

        }
    }
}

TbbSmoothNormalController::TbbSmoothNormalController() {
}

TbbSmoothNormalController::~TbbSmoothNormalController() {
}

void
TbbSmoothNormalController::Synchronize() {
}

}  // end namespace Osd

}  // end namespace OPENSUBDIV_VERSION
}  // end namespace OpenSubdiv
