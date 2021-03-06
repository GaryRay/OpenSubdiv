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

#include <cassert>

#include "../osd/cpuComputeContext.h"
#include "../osd/tbbComputeController.h"
#include "../osd/tbbKernel.h"

#ifdef OPENSUBDIV_HAS_TBB 
    #include <tbb/task_scheduler_init.h>
#endif

namespace OpenSubdiv {
namespace OPENSUBDIV_VERSION {


OsdTbbComputeController::OsdTbbComputeController(int numThreads) {
    _numThreads = numThreads;
    if(_numThreads == -1)
        tbb::task_scheduler_init init;
    else
        tbb::task_scheduler_init init(numThreads);
}


void
OsdTbbComputeController::ApplyBilinearFaceVerticesKernel(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeFace(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::F_IT)->GetBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::F_ITa)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyBilinearEdgeVerticesKernel(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeBilinearEdge(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::E_IT)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyBilinearVertexVerticesKernel(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeBilinearVertex(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyCatmarkFaceVerticesKernel(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeFace(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::F_IT)->GetBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::F_ITa)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyCatmarkEdgeVerticesKernel(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeEdge(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::E_IT)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::E_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyCatmarkVertexVerticesKernelB(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeVertexB(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_IT)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::V_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyCatmarkVertexVerticesKernelA1(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeVertexA(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::V_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd(), false);
}

void
OsdTbbComputeController::ApplyCatmarkVertexVerticesKernelA2(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeVertexA(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::V_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd(), true);
}

void
OsdTbbComputeController::ApplyLoopEdgeVerticesKernel(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeEdge(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::E_IT)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::E_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyLoopVertexVerticesKernelB(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeLoopVertexB(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_IT)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::V_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd());
}

void
OsdTbbComputeController::ApplyLoopVertexVerticesKernelA1(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeVertexA(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::V_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd(), false);
}

void
OsdTbbComputeController::ApplyLoopVertexVerticesKernelA2(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    OsdTbbComputeVertexA(
        context->GetVertexDescriptor(),
        context->GetCurrentVertexBuffer(),
        context->GetCurrentVaryingBuffer(),
        (const int*)context->GetTable(FarSubdivisionTables::V_ITa)->GetBuffer(),
        (const float*)context->GetTable(FarSubdivisionTables::V_W)->GetBuffer(),
        batch.GetVertexOffset(), batch.GetTableOffset(), batch.GetStart(), batch.GetEnd(), true);
}

void
OsdTbbComputeController::ApplyVertexEdits(
    FarKernelBatch const &batch, OsdCpuComputeContext *context) const {

    assert(context);

    const OsdCpuHEditTable *edit = context->GetEditTable(batch.GetTableIndex());
    assert(edit);

    const OsdCpuTable * primvarIndices = edit->GetPrimvarIndices();
    const OsdCpuTable * editValues = edit->GetEditValues();

    if (edit->GetOperation() == FarVertexEdit::Add) {
        OsdTbbEditVertexAdd(context->GetVertexDescriptor(),
                            context->GetCurrentVertexBuffer(),
                            edit->GetPrimvarOffset(),
                            edit->GetPrimvarWidth(),
                            batch.GetVertexOffset(), 
                            batch.GetTableOffset(), 
                            batch.GetStart(), 
                            batch.GetEnd(),
                            static_cast<unsigned int*>(primvarIndices->GetBuffer()),
                            static_cast<float*>(editValues->GetBuffer()));
    } else if (edit->GetOperation() == FarVertexEdit::Set) {
        OsdTbbEditVertexSet(context->GetVertexDescriptor(),
                            context->GetCurrentVertexBuffer(),
                            edit->GetPrimvarOffset(),
                            edit->GetPrimvarWidth(),
                            batch.GetVertexOffset(), 
                            batch.GetTableOffset(), 
                            batch.GetStart(), 
                            batch.GetEnd(),
                            static_cast<unsigned int*>(primvarIndices->GetBuffer()),
                            static_cast<float*>(editValues->GetBuffer()));
    }
}

void
OsdTbbComputeController::Synchronize() {
    // XXX: 
}

}  // end namespace OPENSUBDIV_VERSION
}  // end namespace OpenSubdiv

