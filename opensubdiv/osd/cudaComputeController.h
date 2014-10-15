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

#ifndef OSD_CUDA_COMPUTE_CONTROLLER_H
#define OSD_CUDA_COMPUTE_CONTROLLER_H

#include "../version.h"

#include "../far/kernelBatchDispatcher.h"
#include "../osd/cudaComputeContext.h"
#include "../osd/vertexDescriptor.h"

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {

namespace Osd {

/// \brief Compute controller for launching CUDA subdivision kernels.
///
/// CudaComputeController is a compute controller class to launch
/// Cuda subdivision kernels. It requires CudaVertexBufferInterface
/// as arguments of Refine function.
///
/// Controller entities execute requests from Context instances that they share
/// common interfaces with. Controllers are attached to discrete compute devices
/// and share the devices resources with Context entities.
///
class CudaComputeController {
public:
    typedef CudaComputeContext ComputeContext;

    /// Constructor.
    CudaComputeController();

    /// Destructor.
    ~CudaComputeController();

    /// Execute subdivision kernels and apply to given vertex buffers.
    ///
    /// @param  context       The CudaContext to apply refinement operations to
    ///
    /// @param  batches       Vector of batches of vertices organized by operative
    ///                       kernel
    ///
    /// @param  vertexBuffer  Vertex-interpolated data buffer
    ///
    /// @param  vertexDesc    The descriptor of vertex elements to be refined.
    ///                       if it's null, all primvars in the vertex buffer
    ///                       will be refined.
    ///
    /// @param  varyingBuffer Vertex-interpolated data buffer
    ///
    /// @param  varyingDesc   The descriptor of varying elements to be refined.
    ///                       if it's null, all primvars in the vertex buffer
    ///                       will be refined.
    ///
    template<class VERTEX_BUFFER, class VARYING_BUFFER>
        void Compute( CudaComputeContext const * context,
                      Far::KernelBatchVector const & batches,
                      VERTEX_BUFFER  * vertexBuffer,
                      VARYING_BUFFER * varyingBuffer,
                      VertexBufferDescriptor const * vertexDesc=NULL,
                      VertexBufferDescriptor const * varyingDesc=NULL ){

        if (batches.empty()) return;

        bind(vertexBuffer, varyingBuffer, vertexDesc, varyingDesc);

        Far::KernelBatchDispatcher::Apply(this, context, batches, /*maxlevel*/ -1);

        unbind();
    }

    /// Execute subdivision kernels and apply to given vertex buffers.
    ///
    /// @param  context       The CudaContext to apply refinement operations to
    ///
    /// @param  batches       Vector of batches of vertices organized by operative
    ///                       kernel
    ///
    /// @param  vertexBuffer  Vertex-interpolated data buffer
    ///
    template<class VERTEX_BUFFER>
        void Compute(CudaComputeContext const * context,
                     Far::KernelBatchVector const & batches,
                     VERTEX_BUFFER *vertexBuffer) {

        Compute<VERTEX_BUFFER>(context, batches, vertexBuffer, (VERTEX_BUFFER*)0);
    }

    /// Waits until all running subdivision kernels finish.
    void Synchronize();

protected:

    friend class Far::KernelBatchDispatcher;

    void ApplyStencilTableKernel(Far::KernelBatch const &batch,
        ComputeContext const *context) const;

    template<class VERTEX_BUFFER, class VARYING_BUFFER>
        void bind( VERTEX_BUFFER * vertexBuffer,
                   VARYING_BUFFER * varyingBuffer,
                   VertexBufferDescriptor const * vertexDesc,
                   VertexBufferDescriptor const * varyingDesc ) {

        // if the vertex buffer descriptor is specified, use it.
        // otherwise, assumes the data is tightly packed in the vertex buffer.
        if (vertexDesc) {
            _currentBindState.vertexDesc = *vertexDesc;
        } else {
            int numElements = vertexBuffer ? vertexBuffer->GetNumElements() : 0;
            _currentBindState.vertexDesc =
                VertexBufferDescriptor(0, numElements, numElements);
        }

        if (varyingDesc) {
            _currentBindState.varyingDesc = *varyingDesc;
        } else {
            int numElements = varyingBuffer ? varyingBuffer->GetNumElements() : 0;
            _currentBindState.varyingDesc =
                VertexBufferDescriptor(0, numElements, numElements);
        }

        _currentBindState.vertexBuffer = vertexBuffer ?
            static_cast<float*>(vertexBuffer->BindCudaBuffer()) : 0;
        _currentBindState.varyingBuffer = varyingBuffer ?
            static_cast<float*>(varyingBuffer->BindCudaBuffer()) : 0;
    }

    /// Unbinds any previously bound vertex and varying data buffers.
    void unbind() {
        _currentBindState.Reset();
    }

private:

    // Bind state is a transitional state during refinement.
    // It doesn't take an ownership of the vertex buffers.
    struct BindState {

        BindState() : vertexBuffer(NULL), varyingBuffer(NULL) {}

        void Reset() {
            vertexBuffer = varyingBuffer = NULL;
            vertexDesc.Reset();
            varyingDesc.Reset();
        }

        float *GetVertexBufferAtOffset() const {
            return vertexBuffer ? vertexBuffer + vertexDesc.offset : 0;
        }

        float *GetVaryingBufferAtOffset() const {
            return varyingBuffer ? varyingBuffer + varyingDesc.offset : 0;
        }

        float * vertexBuffer,   // cuda buffers
              * varyingBuffer;

        VertexBufferDescriptor vertexDesc,
                                  varyingDesc;
    };

    BindState _currentBindState;
};

}  // end namespace Osd

}  // end namespace OPENSUBDIV_VERSION
using namespace OPENSUBDIV_VERSION;

}  // end namespace OpenSubdiv

#endif  // OSD_CUDA_COMPUTE_CONTROLLER_H
