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
#include <cmath>
#include <cstring>
#include <cstdio>
#include <cfloat>
#include "convert_bezier.h"

#include "bezier/math.h"
#include "bezier/bezier.h"

using namespace OsdBezier;

static float ef[27] = {
    0.812816f, 0.500000f, 0.363644f, 0.287514f,
    0.238688f, 0.204544f, 0.179229f, 0.159657f,
    0.144042f, 0.131276f, 0.120632f, 0.111614f,
    0.103872f, 0.09715f, 0.0912559f, 0.0860444f,
    0.0814022f, 0.0772401f, 0.0734867f, 0.0700842f,
    0.0669851f, 0.0641504f, 0.0615475f, 0.0591488f,
    0.0569311f, 0.0548745f, 0.0529621f
};

inline float
csf(unsigned int n, unsigned int j)
{
    if (j%2 == 0) {
        return cosf((2.0f * float(M_PI) * float(float(j-0)/2.0f))/(float(n)+3.0f));
    } else {
        return sinf((2.0f * float(M_PI) * float(float(j-1)/2.0f))/(float(n)+3.0f));
    }
}

OsdBezier::matrix4t<float>
computeMatrixSimplified(float sharpness)
{
    float s = pow(2.0f, sharpness);
    float s2 = s*s;
    float s3 = s2*s;

    OsdBezier::matrix4t<float> m = OsdBezier::matrix4t<float>(
        0, s + 1 + 3*s2 - s3, 7*s - 2 - 6*s2 + 2*s3, (1-s)*(s-1)*(s-1),
        0,       (1+s)*(1+s),        6*s - 2 - 2*s2,       (s-1)*(s-1),
        0,               1+s,               6*s - 2,               1-s,
        0,                 1,               6*s - 2,                 1);

    for (int i = 0; i < 16; ++i) {
        m.m[i] /= (s*6.0);
    }
    m[0][0] = 1.0/6.0;

    return m;
}

int convertRegular(std::vector<float> &bezierVertices,
                   std::vector<float> &bezierBounds,
                   std::vector<int> &cpIndices,
                   float const *vertices,
                   OpenSubdiv::Far::PatchTables const *patchTables,
                   OpenSubdiv::Far::PatchTables::PatchArray const &parray)
{
    int numPatches = parray.GetNumPatches();

    // regular to bezier
    for (int i = 0; i < numPatches; i++) {
        const int *verts = &patchTables->GetPatchTable()[parray.GetVertIndex() + i*16];

        // Watertight evaluation.
        // 0   1   2   3
        // 4   5   6   7
        // 8   9  10  11
        // 12  13  14  15
        using namespace OsdBezier;
        vec3f P[16];
        for (int j = 0; j < 16; ++j) {
            P[j] = vec3f(&vertices[verts[j]*3]);
        }
        vec3f BP[16];

        BP[0] = (((P[0] + P[10]) + (P[2] + P[8])) + (4*(P[4] + P[6]) + 4*(P[1] + P[9])))+ 16*P[5];
        BP[1] = (4*(P[4] + P[6]) + 2*(P[8] + P[10])) + 4*(4*P[5] + 2*P[9]);
        BP[2] = (4*(P[8] + P[10]) + 2*(P[4] + P[6])) + 4*(4*P[9] + 2*P[5]);
        BP[3] = (((P[12] + P[6]) + (P[4] + P[14])) + (4*(P[8] + P[10]) + 4*(P[5] + P[13]))) + 16*P[9];
        BP[4] = (4*(P[1] + P[9]) + 2*(P[2] + P[10])) + 4*(4*P[5] + 2*P[6]);
        BP[5] = (4*(4*P[5] + 2*P[6]) + 2*(2*P[10] + 4*P[9]));
        BP[6] = (2*(4*P[5] + 2*P[6]) + 4*(2*P[10] + 4*P[9]));
        BP[7] = (4*(P[5] + P[13]) + 2*(P[6] + P[14])) + 4*(4*P[9] + 2*P[10]);
        BP[8] = (4*(P[2] + P[10]) + 2*(P[1] + P[9])) + 4*(2*P[5] + 4*P[6]);
        BP[9] = (4*(2*P[5] + 4*P[6]) + 2*(4*P[10] + 2*P[9]));
        BP[10] = (2*(2*P[5] + 4*P[6]) + 4*(4*P[10] + 2*P[9]));
        BP[11] = (4*(P[6] + P[14]) + 2*(P[5] + P[13])) + 4*(4*P[10] + 2*P[9]);
        BP[12] = (((P[3] + P[9]) + (P[1] + P[11])) + (4*(P[2] + P[10]) + 4*(P[5] + P[7]))) + 16*P[6];
        BP[13] = (4*(P[5] + P[7]) + 2*(P[9] + P[11])) + 4*(4*P[6] + 2*P[10]);
        BP[14] = (4*(P[9] + P[11]) + 2*(P[5] + P[7])) + 4*(4*P[10] + 2*P[6]);
        BP[15] =(((P[5] + P[15]) + (P[7] + P[13])) + (4*(P[6]+P[14]) + 4*(P[9]+P[11]))) + 16*P[10];
        for (int j = 0; j < 16; ++j) {
            BP[j] = BP[j]*(1.0f/36.0f);
        }
        vec3f min(FLT_MAX), max(-FLT_MAX);
        for (int j = 0; j < 16; ++j) {
            min = min.min(BP[j]);
            max = max.max(BP[j]);
            bezierVertices.push_back(BP[j][0]);
            bezierVertices.push_back(BP[j][1]);
            bezierVertices.push_back(BP[j][2]);
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);

        // save center quad indices
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 5]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 6]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 9]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 10]);
    }
    return numPatches;
}

int convertSingleCrease(std::vector<float> &bezierVertices,
                        std::vector<float> &bezierBounds,
                        std::vector<int> &cpIndices,
                        float const *vertices,
                        OpenSubdiv::Far::PatchTables const *patchTables,
                        OpenSubdiv::Far::PatchTables::PatchArray const &parray)
{
    int numPatches = parray.GetNumPatches();

    // single-crease to bezier.
    for (int i = 0; i < numPatches; i++) {
        const int *verts = &patchTables->GetPatchTable()[parray.GetVertIndex() + i*16];

        // same watertight bezier conversion should be applied as regular patch.

        // Watertight Splinet to Bezier evaluation.
        using namespace OsdBezier;
        vec3f P[16];
        for (int j = 0; j < 16; ++j) {
            P[j] = vec3f(&vertices[verts[j]*3]);
        }
        // sharpness matrix
        int index = patchTables->GetSharpnessIndexTable()[parray.GetPatchIndex() + i];
        float sharpness = patchTables->GetSharpnessValues()[index];
        //float Sf = floor(sharpness);
        //float Sc = ceil(sharpness);
        //float Sr = sharpness - Sf;
        //matrix4t<float> Mf = computeMatrixSimplified(Sf);
        //matrix4t<float> Mc = computeMatrixSimplified(Sc);
        //matrix4t<float> Mj = (1-Sr) * Mf + Sr * Mi;

        matrix4t<float> Ms = computeMatrixSimplified(sharpness);

        // convert bspline to bezier
        const float Q[4][4] = {
            { 1.f/6.f, 4.f/6.f, 1.f/6.f, 0.f },
            { 0.f,     4.f/6.f, 2.f/6.f, 0.f },
            { 0.f,     2.f/6.f, 4.f/6.f, 0.f },
            { 0.f,     1.f/6.f, 4.f/6.f, 1.f/6.f } };

        // Infinite sharp
        const float Mi[4][4] = {
            { 1.f/6.f, 4.f/6.f, 1.f/6.f, 0.f },
            { 0.f,     4.f/6.f, 2.f/6.f, 0.f },
            { 0.f,     2.f/6.f, 4.f/6.f, 0.f },
            { 0.f,     0.f,     1.f,     0.f } };

        // ----------------------------------------------------
        vec3f cp0[16], cp1[16];
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                vec3f H[4];
                for (int l = 0; l < 4; l++) {
                    H[l] = vec3f(0.0f);
                    for (int m = 0; m < 4; m++) {
                        int vert = verts[l*4 + m];
                        H[l] += Q[j][m] * vec3f(vertices[vert*3+0],
                                                vertices[vert*3+1],
                                                vertices[vert*3+2]);
                    }
                }
                for (int m = 0; m < 4; m++) {
                    cp0[j*4+k] += Ms[k][m] * H[m];
                    cp1[j*4+k] += Mi[k][m] * H[m];
                }
            }
        }
        // split patch
        BezierPatch<vec3f, float, 4> patch0(cp0);
        BezierPatch<vec3f, float, 4> patch1(cp1);

        BezierPatch<vec3f, float, 4> patch0s[2];
        BezierPatch<vec3f, float, 4> patch1s[2];
        float t = 1.0 - pow(2.0f, -sharpness);
        patch0.SplitU(patch0s, t);
        patch1.SplitU(patch1s, t);

        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                cp0[j*4+k] = patch0s[1].Get(k, j);
                cp1[j*4+k] = patch1s[0].Get(k, j);
            }
        }

        vec3f min(FLT_MAX), max(-FLT_MAX);
        for (int j = 0; j < 16; ++j) {
            bezierVertices.push_back(cp0[j][0]);
            bezierVertices.push_back(cp0[j][1]);
            bezierVertices.push_back(cp0[j][2]);
            min = min.min(cp0[j]);
            max = max.max(cp0[j]);
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);
        min = vec3f(FLT_MAX);
        max = vec3f(-FLT_MAX);

        for (int j = 0; j < 16; ++j) {
            bezierVertices.push_back(cp1[j][0]);
            bezierVertices.push_back(cp1[j][1]);
            bezierVertices.push_back(cp1[j][2]);
            min = min.min(cp1[j]);
            max = max.max(cp1[j]);
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);

#if 0
        // save center quad indices
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 5]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 6]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 9]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 10]);
        // save center quad indices
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 5]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 6]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 9]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*16 + 10]);
#endif
        for (int j = 0; j < 8; ++j) {
            cpIndices.push_back(-1);
        }

    }
    return numPatches*2;
}

int convertBoundary(std::vector<float> &bezierVertices,
                    std::vector<float> &bezierBounds,
                    std::vector<int> &cpIndices,
                    float const *vertices,
                    OpenSubdiv::Far::PatchTables const *patchTables,
                    OpenSubdiv::Far::PatchTables::PatchArray const &parray)
{
    int numPatches = parray.GetNumPatches();

    // convert boundary bspline to bezier
    const float Q[4][4] = {
        { 1.f/6.f, 4.f/6.f, 1.f/6.f, 0.f },
        { 0.f,     4.f/6.f, 2.f/6.f, 0.f },
        { 0.f,     2.f/6.f, 4.f/6.f, 0.f },
        { 0.f,     1.f/6.f, 4.f/6.f, 1.f/6.f } };
    const float B[4][3] = {
        {1.f,     0.f,     0.f},
        {4.f/6.f, 2.f/6.f, 0.f},
        {2.f/6.f, 4.f/6.f, 0.f},
        {1.f/6.f, 4.f/6.f, 1.f/6.f}};

    // regular to bezier
    for (int i = 0; i < numPatches; i++) {
        vec3f min(FLT_MAX), max(-FLT_MAX);
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                vec3f H[3];
                for (int l = 0; l < 3; l++) {
                    H[l] = vec3f(0.0f);
                    for (int m = 0; m < 4; m++) {
                        int vert = patchTables->GetPatchTable()[parray.GetVertIndex() + i*12 + l*4 + m];
                        H[l] += Q[j][m] * vec3f(vertices[vert*3+0],
                                                vertices[vert*3+1],
                                                vertices[vert*3+2]);
                    }
                }
                vec3f cp(0.0f);
                for (int m = 0; m < 3; m++) {
                    cp += B[k][m] * H[m];
                }
                min = min.min(cp);
                max = max.max(cp);
                bezierVertices.push_back(cp[0]);
                bezierVertices.push_back(cp[1]);
                bezierVertices.push_back(cp[2]);
            }
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);

        // save center quad indices
#if 1
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*12 + 1]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*12 + 2]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*12 + 5]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*12 + 6]);
#else
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
#endif
    }
    return numPatches;
}

int convertCorner(std::vector<float> &bezierVertices,
                  std::vector<float> &bezierBounds,
                  std::vector<int> &cpIndices,
                  float const *vertices,
                  OpenSubdiv::Far::PatchTables const *patchTables,
                  OpenSubdiv::Far::PatchTables::PatchArray const &parray)
{
    int numPatches = parray.GetNumPatches();

    // convert bspline to bezier
    const float B[4][3] = {
        {1.f,     0.f,     0.f},
        {4.f/6.f, 2.f/6.f, 0.f},
        {2.f/6.f, 4.f/6.f, 0.f},
        {1.f/6.f, 4.f/6.f, 1.f/6.f}};

    // regular to bezier
    for (int i = 0; i < numPatches; i++) {
        vec3f min(FLT_MAX), max(-FLT_MAX);
        for (int j = 0; j < 4; j++) {
            for (int k = 0; k < 4; k++) {
                float H[3][3];
                for (int l = 0; l < 3; l++) {
                    H[l][0] = H[l][1] = H[l][2] = 0;
                    for (int m = 0; m < 3; m++) {
                        int vert = patchTables->GetPatchTable()[parray.GetVertIndex() + i*9 + l*3 + m];
                        H[l][0] += B[3-j][2-m] * vertices[vert*3+0];
                        H[l][1] += B[3-j][2-m] * vertices[vert*3+1];
                        H[l][2] += B[3-j][2-m] * vertices[vert*3+2];
                    }
                }
                float cp[3];
                cp[0] = cp[1] = cp[2] = 0;
                for (int m = 0; m < 3; m++) {
                    cp[0] += B[k][m] * H[m][0];
                    cp[1] += B[k][m] * H[m][1];
                    cp[2] += B[k][m] * H[m][2];
                }
                min = min.min(cp);
                max = max.max(cp);
                bezierVertices.push_back(cp[0]);
                bezierVertices.push_back(cp[1]);
                bezierVertices.push_back(cp[2]);
            }
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);

        // save center quad indices
#if 0
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*9 + 0]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*9 + 3]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*9 + 1]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + i*9 + 4]);
#else
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
#endif
    }
    return numPatches;
}


int convertGregory(std::vector<float> &bezierVertices,
                   std::vector<float> &bezierBounds,
                   std::vector<int> &cpIndices,
                   float const *vertices,
                   OpenSubdiv::Far::PatchTables const *patchTables,
                   OpenSubdiv::Far::PatchTables::PatchArray const &parray)
{
    int const * vertexValenceBuffer = &patchTables->GetVertexValenceTable()[0];
    int valences[4];
    int length = 3;
    int stride = 3;
    int maxValence = patchTables->GetMaxValence();
    printf("Max Valence = %d\n", maxValence);
    float *r = (float*)alloca((maxValence+2)*4*length*sizeof(float));
    float *f = (float*)alloca(maxValence*length*sizeof(float)),
        *pos = (float*)alloca(length*sizeof(float)),
        *opos = (float*)alloca(length*4*sizeof(float));
    float const * inOffset = (float const*)&vertices[0];

    float *Ep=(float*)alloca(length*4*sizeof(float)),
        *Em=(float*)alloca(length*4*sizeof(float)),
        *Fp=(float*)alloca(length*4*sizeof(float)),
        *Fm=(float*)alloca(length*4*sizeof(float));
    float *Em_ip=(float*)alloca(length*sizeof(float)),
        *Ep_im=(float*)alloca(length*sizeof(float));
    float *q=(float*)alloca(length*16*sizeof(float));

    for (int patchIndex = 0; patchIndex < (int)parray.GetNumPatches(); ++patchIndex) {
        float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
        float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

        int const * vertexIndices = &patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4];
        int const *quadOffsetBuffer = &patchTables->GetQuadOffsetTable()[parray.GetQuadOffsetIndex() + patchIndex*4];

        bool badPatch = false;

        float *rp,
            *e0 = r + maxValence*4*length,
            *e1 = e0 + 4*length;
        memset(r, 0, (maxValence+2)*4*length*sizeof(float));

        memset(opos, 0, length*4*sizeof(float));

        for (int vid=0; vid < 4; ++vid) {

            int vertexID = vertexIndices[vid];

            const int *valenceTable = vertexValenceBuffer + vertexID * (2*maxValence+1);
            int valence = abs(*valenceTable);
            assert(valence<=maxValence);
            valences[vid] = valence;

            memcpy(pos, inOffset + vertexID*stride, length*sizeof(float));

            rp=r+vid*maxValence*length;

            int vofs = vid*length;

            for (int i=0; i<valence; ++i) {
                unsigned int im = (i+valence-1)%valence,
                    ip = (i+1)%valence;

                int idx_neighbor   = valenceTable[2*i  + 0 + 1];
                int idx_diagonal   = valenceTable[2*i  + 1 + 1];
                int idx_neighbor_p = valenceTable[2*ip + 0 + 1];
                int idx_neighbor_m = valenceTable[2*im + 0 + 1];
                int idx_diagonal_m = valenceTable[2*im + 1 + 1];

                float const * neighbor   = inOffset + idx_neighbor   * stride;
                float const * diagonal   = inOffset + idx_diagonal   * stride;
                float const * neighbor_p = inOffset + idx_neighbor_p * stride;
                float const * neighbor_m = inOffset + idx_neighbor_m * stride;
                float const * diagonal_m = inOffset + idx_diagonal_m * stride;

                float  *fp = f+i*length;

                for (int k=0; k<length; ++k) {
                    fp[k] = (pos[k]*float(valence) + (neighbor_p[k]+neighbor[k])*2.0f + diagonal[k])/(float(valence)+5.0f);

                    opos[vofs+k] += fp[k];
                    rp[i*length+k] =(neighbor_p[k]-neighbor_m[k])/3.0f + (diagonal[k]-diagonal_m[k])/6.0f;
                }

            }

            for (int k=0; k<length; ++k) {
                opos[vofs+k] /= valence;
            }

            if (valence >= 3 and valence < 30) {
                for (int i=0; i<valence; ++i) {
                    int im = (i+valence-1)%valence;
                    for (int k=0; k<length; ++k) {
                        float e = 0.5f*(f[i*length+k]+f[im*length+k]);
                        e0[vofs+k] += csf(valence-3, 2*i) * e;
                        e1[vofs+k] += csf(valence-3, 2*i+1) * e;
                    }
                }
                for (int k=0; k<length; ++k) {
                    e0[vofs+k] *= ef[valence-3];
                    e1[vofs+k] *= ef[valence-3];
                }
            } else {
                badPatch = true;
            }
        }

        // tess control

        // Control Vertices based on :
        // "Approximating Subdivision Surfaces with Gregory Patches for Hardware Tessellation"
        // Loop, Schaefer, Ni, Castafio (ACM ToG Siggraph Asia 2009)
        //
        //  P3         e3-      e2+         E2
        //     O--------O--------O--------O
        //     |        |        |        |
        //     |        |        |        |
        //     |        | f3-    | f2+    |
        //     |        O        O        |
        // e3+ O------O            O------O e2-
        //     |     f3+          f2-     |
        //     |                          |
        //     |                          |
        //     |      f0-         f1+     |
        // e0- O------O            O------O e1+
        //     |        O        O        |
        //     |        | f0+    | f1-    |
        //     |        |        |        |
        //     |        |        |        |
        //     O--------O--------O--------O
        //  P0         e0+      e1-         E1
        //

        for (int vid=0; vid<4; ++vid) {

            int ip = (vid+1)%4;
            int im = (vid+3)%4;
            int n = valences[vid];
            int const *quadOffsets = quadOffsetBuffer;

            int start = quadOffsets[vid] & 0x00ff;
            int prev = (quadOffsets[vid] & 0xff00) / 256;

            for (int k=0, ofs=vid*length; k<length; ++k, ++ofs) {

                Ep[ofs] = opos[ofs] + e0[ofs] * csf(n-3, 2*start) + e1[ofs]*csf(n-3, 2*start +1);
                Em[ofs] = opos[ofs] + e0[ofs] * csf(n-3, 2*prev ) + e1[ofs]*csf(n-3, 2*prev + 1);

            }

            unsigned int np = valences[ip],
                nm = valences[im];

            unsigned int prev_p = (quadOffsets[ip] & 0xff00) / 256,
                start_m = quadOffsets[im] & 0x00ff;

            for (int k=0, ipofs=ip*length, imofs=im*length; k<length; ++k, ++ipofs, ++imofs) {
                Em_ip[k] = opos[ipofs] + e0[ipofs]*csf(np-3, 2*prev_p)  + e1[ipofs]*csf(np-3, 2*prev_p+1);
                Ep_im[k] = opos[imofs] + e0[imofs]*csf(nm-3, 2*start_m) + e1[imofs]*csf(nm-3, 2*start_m+1);
            }

            float s1 = 3.0f - 2.0f*csf(n-3,2)-csf(np-3,2),
                s2 = 2.0f*csf(n-3,2),
                s3 = 3.0f -2.0f*cosf(2.0f*float(M_PI)/float(n)) - cosf(2.0f*float(M_PI)/float(nm));

            rp = r + vid*maxValence*length;
            for (int k=0, ofs=vid*length; k<length; ++k, ++ofs) {
                Fp[ofs] = (csf(np-3,2)*opos[ofs] + s1*Ep[ofs] + s2*Em_ip[k] + rp[start*length+k])/3.0f;
                Fm[ofs] = (csf(nm-3,2)*opos[ofs] + s3*Em[ofs] + s2*Ep_im[k] - rp[prev*length+k])/3.0f;
            }
        }

        float * p[20];
        for (int i=0, ofs=0; i<4; ++i, ofs+=length) {
            p[i*5+0] = opos + ofs;
            p[i*5+1] =   Ep + ofs;
            p[i*5+2] =   Em + ofs;
            p[i*5+3] =   Fp + ofs;
            p[i*5+4] =   Fm + ofs;
        }

        float u = 0.5f, v=0.5f;
        float U = 1-u, V=1-v;
        float d11 = u+v; if(u+v==0.0f) d11 = 1.0f;
        float d12 = U+v; if(U+v==0.0f) d12 = 1.0f;
        float d21 = u+V; if(u+V==0.0f) d21 = 1.0f;
        float d22 = U+V; if(U+V==0.0f) d22 = 1.0f;

        for (int k=0; k<length; ++k) {
            q[ 5*length+k] = (u*p[ 3][k] + v*p[ 4][k])/d11;
            q[ 6*length+k] = (U*p[ 9][k] + v*p[ 8][k])/d12;
            q[ 9*length+k] = (u*p[19][k] + V*p[18][k])/d21;
            q[10*length+k] = (U*p[13][k] + V*p[14][k])/d22;
        }

        memcpy(q+ 0*length, p[ 0], length*sizeof(float));
        memcpy(q+ 1*length, p[ 1], length*sizeof(float));
        memcpy(q+ 2*length, p[ 7], length*sizeof(float));
        memcpy(q+ 3*length, p[ 5], length*sizeof(float));
        memcpy(q+ 4*length, p[ 2], length*sizeof(float));
        memcpy(q+ 7*length, p[ 6], length*sizeof(float));
        memcpy(q+ 8*length, p[16], length*sizeof(float));
        memcpy(q+11*length, p[12], length*sizeof(float));
        memcpy(q+12*length, p[15], length*sizeof(float));
        memcpy(q+13*length, p[17], length*sizeof(float));
        memcpy(q+14*length, p[11], length*sizeof(float));
        memcpy(q+15*length, p[10], length*sizeof(float));

        if (badPatch) {
            memset(q, 0, 16*3*sizeof(float));
        }

        for (int x = 0; x < 4; ++x) {
            for (int y = 0; y < 4; ++y) {
                int index = y*4 + x;
                min[0] = std::min(min[0], q[index*3+0]);
                min[1] = std::min(min[1], q[index*3+1]);
                min[2] = std::min(min[2], q[index*3+2]);
                max[0] = std::max(max[0], q[index*3+0]);
                max[1] = std::max(max[1], q[index*3+1]);
                max[2] = std::max(max[2], q[index*3+2]);
                bezierVertices.push_back(q[index*3+0]);
                bezierVertices.push_back(q[index*3+1]);
                bezierVertices.push_back(q[index*3+2]);
            }
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);


        // save center quad indices
        if (!badPatch) {
            cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 0]);
            cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 1]);
            cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 2]);
            cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 3]);
        } else {
            cpIndices.push_back(-1);
            cpIndices.push_back(-1);
            cpIndices.push_back(-1);
            cpIndices.push_back(-1);
        }
    }

    return parray.GetNumPatches();
}

int convertBoundaryGregory(std::vector<float> &bezierVertices,
                           std::vector<float> &bezierBounds,
                           std::vector<int> &cpIndices,
                           float const *vertices,
                           OpenSubdiv::Far::PatchTables const *patchTables,
                           OpenSubdiv::Far::PatchTables::PatchArray const &parray)
{
    int const * vertexValenceBuffer = &patchTables->GetVertexValenceTable()[0];
    int valences[4], zerothNeighbors[4];
    int length = 3;
    int stride = 3;
    int maxValence = patchTables->GetMaxValence();
    float const * inOffset = (float const*)&vertices[0];
    float  *r  = (float*)alloca((maxValence+2)*4*length*sizeof(float));
    float *f=(float*)alloca(maxValence*length*sizeof(float)),
        *org=(float*)alloca(length*4*sizeof(float)),
        *opos=(float*)alloca(length*4*sizeof(float));
    float *Ep=(float*)alloca(length*4*sizeof(float)),
        *Em=(float*)alloca(length*4*sizeof(float)),
        *Fp=(float*)alloca(length*4*sizeof(float)),
        *Fm=(float*)alloca(length*4*sizeof(float));
    float *Em_ip=(float*)alloca(length*sizeof(float)),
        *Ep_im=(float*)alloca(length*sizeof(float));
    float *q=(float*)alloca(length*16*sizeof(float));

    for (int patchIndex = 0; patchIndex < (int)parray.GetNumPatches(); ++patchIndex) {
        float min[3] = {FLT_MAX, FLT_MAX, FLT_MAX};
        float max[3] = {-FLT_MAX, -FLT_MAX, -FLT_MAX};

        // vertex
        int const * vertexIndices = &patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4];
        int const * quadOffsetBuffer = &patchTables->GetQuadOffsetTable()[parray.GetQuadOffsetIndex() + patchIndex*4];

        bool badPatch = false;

        float *rp,
            *e0 = r + maxValence*4*length,
            *e1 = e0 + 4*length;
        memset(r, 0, (maxValence+2)*4*length*sizeof(float));

        memset(opos, 0, length*4*sizeof(float));

        for (int vid=0; vid < 4; ++vid) {

            int vertexID = vertexIndices[vid];

            const int *valenceTable = vertexValenceBuffer + vertexID * (2*maxValence+1);
            int valence = *valenceTable,
                ivalence = abs(valence);

            assert(ivalence<=maxValence);
            valences[vid] = valence;

            int vofs = vid * length;

            float *pos=org + vofs;
            memcpy(pos, inOffset + vertexID*stride, length*sizeof(float));

            int boundaryEdgeNeighbors[2];
            unsigned int currNeighbor = 0,
                ibefore=0,
                zerothNeighbor=0;

            rp=r+vid*maxValence*length;

            for (int i=0; i<ivalence; ++i) {
                unsigned int im = (i+ivalence-1)%ivalence,
                    ip = (i+1)%ivalence;

                int idx_neighbor   = valenceTable[2*i  + 0 + 1];
                int idx_diagonal   = valenceTable[2*i  + 1 + 1];
                int idx_neighbor_p = valenceTable[2*ip + 0 + 1];
                int idx_neighbor_m = valenceTable[2*im + 0 + 1];
                int idx_diagonal_m = valenceTable[2*im + 1 + 1];

                int valenceNeighbor = vertexValenceBuffer[idx_neighbor * (2*maxValence+1)];
                if (valenceNeighbor < 0) {

                    if (currNeighbor<2) {
                        boundaryEdgeNeighbors[currNeighbor] = idx_neighbor;
                    }
                    currNeighbor++;

                    if (currNeighbor == 1)    {
                        ibefore = i;
                        zerothNeighbor = i;
                    } else {
                        if (i-ibefore == 1) {
                            int tmp = boundaryEdgeNeighbors[0];
                            boundaryEdgeNeighbors[0] = boundaryEdgeNeighbors[1];
                            boundaryEdgeNeighbors[1] = tmp;
                            zerothNeighbor = i;
                        }
                    }
                }

                float const * neighbor   = inOffset + idx_neighbor   * stride;
                float const * diagonal   = inOffset + idx_diagonal   * stride;
                float const * neighbor_p = inOffset + idx_neighbor_p * stride;
                float const * neighbor_m = inOffset + idx_neighbor_m * stride;
                float const * diagonal_m = inOffset + idx_diagonal_m * stride;

                float *fp = f+i*length;

                for (int k=0; k<length; ++k) {
                    fp[k] = (pos[k]*float(ivalence) + (neighbor_p[k]+neighbor[k])*2.0f + diagonal[k])/(float(ivalence)+5.0f);

                    opos[vofs+k] += fp[k];
                    rp[i*length+k] =(neighbor_p[k]-neighbor_m[k])/3.0f + (diagonal[k]-diagonal_m[k])/6.0f;
                }
            }

            for (int k=0; k<length; ++k) {
                opos[vofs+k] /= ivalence;
            }

            zerothNeighbors[vid] = zerothNeighbor;

            if (currNeighbor == 1) {
                boundaryEdgeNeighbors[1] = boundaryEdgeNeighbors[0];
            }

            if (ivalence >= 3) {
                for (int i=0; i<ivalence; ++i) {
                    unsigned int im = (i+ivalence-1)%ivalence;
                    for (int k=0; k<length; ++k) {
                        float e = 0.5f*(f[i*length+k]+f[im*length+k]);
                        e0[vofs+k] += csf(ivalence-3, 2*i  ) * e;
                        e1[vofs+k] += csf(ivalence-3, 2*i+1) * e;
                    }
                }
                for (int k=0; k<length; ++k) {
                    e0[vofs+k] *= ef[ivalence-3];
                    e1[vofs+k] *= ef[ivalence-3];
                }
            } else {
                badPatch = true;
            }

            if (valence<0) {
                if (ivalence>2) {
                    for (int k=0; k<length; ++k) {
                        opos[vofs+k] = (inOffset[boundaryEdgeNeighbors[0]*stride+k] +
                                        inOffset[boundaryEdgeNeighbors[1]*stride+k] + 4.0f*pos[k])/6.0f;
                    }
                } else {
                    memcpy(opos, pos, length*sizeof(float));
                }

                float k = float(float(ivalence) - 1.0f);    //k is the number of faces
                float c = cosf(float(M_PI)/k);
                float s = sinf(float(M_PI)/k);
                float gamma = -(4.0f*s)/(3.0f*k+c);
                float alpha_0k = -((1.0f+2.0f*c)*sqrtf(1.0f+c))/((3.0f*k+c)*sqrtf(1.0f-c));
                float beta_0 = s/(3.0f*k + c);

                int idx_diagonal = valenceTable[2*zerothNeighbor + 1 + 1];
                assert(idx_diagonal>0);
                float const * diagonal = inOffset + idx_diagonal * stride;

                for (int j=0; j<length; ++j) {
                    e0[vofs+j] = (inOffset[boundaryEdgeNeighbors[0]*stride+j] -
                                  inOffset[boundaryEdgeNeighbors[1]*stride+j])/6.0f;

                    e1[vofs+j] = gamma * pos[j] + beta_0 * diagonal[j] +
                        (inOffset[boundaryEdgeNeighbors[0]*stride+j] +
                         inOffset[boundaryEdgeNeighbors[1]*stride+j]) * alpha_0k;

                }

                for (int x=1; x<ivalence-1; ++x) {
                    unsigned int curri = ((x + zerothNeighbor)%ivalence);
                    float alpha = (4.0f*sinf((float(M_PI) * float(x))/k))/(3.0f*k+c);
                    float beta = (sinf((float(M_PI) * float(x))/k) + sinf((float(M_PI) * float(x+1))/k))/(3.0f*k+c);

                    int idx_neighbor = valenceTable[2*curri + 0 + 1];
                    idx_diagonal = valenceTable[2*curri + 1 + 1];
                    assert( idx_neighbor>0 and idx_diagonal>0 );

                    float const * neighbor = inOffset + idx_neighbor * stride;
                    diagonal = inOffset + idx_diagonal * stride;

                    for (int j=0; j<length; ++j) {
                        e1[vofs+j] += alpha*neighbor[j] + beta*diagonal[j];
                    }
                }

                for (int j=0; j<length; ++j) {
                    e1[vofs+j] /= 3.0f;
                }
            }
        }

        // tess control

        // Control Vertices based on :
        // "Approximating Subdivision Surfaces with Gregory Patches for Hardware Tessellation"
        // Loop, Schaefer, Ni, Castafio (ACM ToG Siggraph Asia 2009)
        //
        //  P3         e3-      e2+         E2
        //     O--------O--------O--------O
        //     |        |        |        |
        //     |        |        |        |
        //     |        | f3-    | f2+    |
        //     |        O        O        |
        // e3+ O------O            O------O e2-
        //     |     f3+          f2-     |
        //     |                          |
        //     |                          |
        //     |      f0-         f1+     |
        // e0- O------O            O------O e1+
        //     |        O        O        |
        //     |        | f0+    | f1-    |
        //     |        |        |        |
        //     |        |        |        |
        //     O--------O--------O--------O
        //  P0         e0+      e1-         E1
        //

        for (int vid=0; vid<4; ++vid) {

            unsigned int ip = (vid+1)%4,
                im = (vid+3)%4,
                n = abs(valences[vid]),
                ivalence = n;

            const int *quadOffsets = quadOffsetBuffer;

            int vofs = vid * length;

            unsigned int   start =  quadOffsets[vid] & 0x00ff,
                prev = (quadOffsets[vid] & 0xff00) / 256,
                np = abs(valences[ip]),
                nm = abs(valences[im]),
                start_m =  quadOffsets[im] & 0x00ff,
                prev_p = (quadOffsets[ip] & 0xff00) / 256;

            if (valences[ip]<-2) {
                unsigned int j = (np + prev_p - zerothNeighbors[ip]) % np;
                for (int k=0, ipofs=ip*length; k<length; ++k, ++ipofs) {
                    Em_ip[k] = opos[ipofs] + cosf((float(M_PI)*j)/float(np-1))*e0[ipofs] + sinf((float(M_PI)*j)/float(np-1))*e1[ipofs];
                }
            } else {
                for (int k=0, ipofs=ip*length; k<length; ++k, ++ipofs) {
                    Em_ip[k] = opos[ipofs] + e0[ipofs]*csf(np-3,2*prev_p)  + e1[ipofs]*csf(np-3,2*prev_p+1);
                }
            }

            if (valences[im]<-2) {
                unsigned int j = (nm + start_m - zerothNeighbors[im]) % nm;
                for (int k=0, imofs=im*length; k<length; ++k, ++imofs) {
                    Ep_im[k] = opos[imofs] + cosf((float(M_PI)*j)/float(nm-1))*e0[imofs] + sinf((float(M_PI)*j)/float(nm-1))*e1[imofs];
                }
            } else {
                for (int k=0, imofs=im*length; k<length; ++k, ++imofs) {
                    Ep_im[k] = opos[imofs] + e0[imofs]*csf(nm-3,2*start_m) + e1[imofs]*csf(nm-3,2*start_m+1);
                }
            }

            if (valences[vid] < 0) {
                n = (n-1)*2;
            }
            if (valences[im] < 0) {
                nm = (nm-1)*2;
            }
            if (valences[ip] < 0) {
                np = (np-1)*2;
            }

            rp=r+vid*maxValence*length;

            if (valences[vid] > 2) {
                float s1 = 3.0f - 2.0f*csf(n-3,2)-csf(np-3,2),
                    s2 = 2.0f*csf(n-3,2),
                    s3 = 3.0f -2.0f*cosf(2.0f*float(M_PI)/float(n)) - cosf(2.0f*float(M_PI)/float(nm));

                for (int k=0, ofs=vofs; k<length; ++k, ++ofs) {
                    Ep[ofs] = opos[ofs] + e0[ofs] * csf(n-3, 2*start) + e1[ofs]*csf(n-3, 2*start +1);
                    Em[ofs] = opos[ofs] + e0[ofs] * csf(n-3, 2*prev ) + e1[ofs]*csf(n-3, 2*prev + 1);
                    Fp[ofs] = (csf(np-3,2)*opos[ofs] + s1*Ep[ofs] + s2*Em_ip[k] + rp[start*length+k])/3.0f;
                    Fm[ofs] = (csf(nm-3,2)*opos[ofs] + s3*Em[ofs] + s2*Ep_im[k] - rp[prev*length+k])/3.0f;
                }
            } else if (valences[vid] < -2) {
                unsigned int jp = (ivalence + start - zerothNeighbors[vid]) % ivalence,
                    jm = (ivalence + prev  - zerothNeighbors[vid]) % ivalence;

                float s1 = 3-2*csf(n-3,2)-csf(np-3,2),
                    s2 = 2*csf(n-3,2),
                    s3 = 3.0f-2.0f*cosf(2.0f*float(M_PI)/n)-cosf(2.0f*float(M_PI)/nm);

                for (int k=0, ofs=vofs; k<length; ++k, ++ofs) {
                    Ep[ofs] = opos[ofs] + cosf((float(M_PI)*jp)/float(ivalence-1))*e0[ofs] + sinf((float(M_PI)*jp)/float(ivalence-1))*e1[ofs];
                    Em[ofs] = opos[ofs] + cosf((float(M_PI)*jm)/float(ivalence-1))*e0[ofs] + sinf((float(M_PI)*jm)/float(ivalence-1))*e1[ofs];
                    Fp[ofs] = (csf(np-3,2)*opos[ofs] + s1*Ep[ofs] + s2*Em_ip[k] + rp[start*length+k])/3.0f;
                    Fm[ofs] = (csf(nm-3,2)*opos[ofs] + s3*Em[ofs] + s2*Ep_im[k] - rp[prev*length+k])/3.0f;
                }

                if (valences[im]<0) {
                    s1=3-2*csf(n-3,2)-csf(np-3,2);
                    for (int k=0, ofs=vofs; k<length; ++k, ++ofs) {
                        Fp[ofs] = Fm[ofs] = (csf(np-3,2)*opos[ofs] + s1*Ep[ofs] + s2*Em_ip[k] + rp[start*length+k])/3.0f;
                    }
                } else if (valences[ip]<0) {
                    s1 = 3.0f-2.0f*cosf(2.0f*float(M_PI)/n)-cosf(2.0f*float(M_PI)/nm);
                    for (int k=0, ofs=vofs; k<length; ++k, ++ofs) {
                        Fm[ofs] = Fp[ofs] = (csf(nm-3,2)*opos[ofs] + s1*Em[ofs] + s2*Ep_im[k] - rp[prev*length+k])/3.0f;
                    }
                }
            } else if (valences[vid]==-2) {
                for (int k=0, ofs=vofs, ipofs=ip*length, imofs=im*length; k<length; ++k, ++ofs, ++ipofs, ++imofs) {
                    Ep[ofs] = (2.0f * org[ofs] + org[ipofs])/3.0f;
                    Em[ofs] = (2.0f * org[ofs] + org[imofs])/3.0f;
                    Fp[ofs] = Fm[ofs] = (4.0f * org[ofs] + org[((vid+2)%n)*stride+k] + 2.0f * org[ipofs] + 2.0f * org[imofs])/9.0f;
                }
            }
        }

        float * p[20];
        for (int vid=0, ofs=0; vid<4; ++vid, ofs+=length) {
            p[vid*5+0] = opos + ofs;
            p[vid*5+1] =   Ep + ofs;
            p[vid*5+2] =   Em + ofs;
            p[vid*5+3] =   Fp + ofs;
            p[vid*5+4] =   Fm + ofs;
        }

        float u = 0.5f, v = 0.5f;
        float U = 1-u, V=1-v;
        float d11 = u+v; if(u+v==0.0f) d11 = 1.0f;
        float d12 = U+v; if(U+v==0.0f) d12 = 1.0f;
        float d21 = u+V; if(u+V==0.0f) d21 = 1.0f;
        float d22 = U+V; if(U+V==0.0f) d22 = 1.0f;

        for (int k=0; k<length; ++k) {
            q[ 5*length+k] = (u*p[ 3][k] + v*p[ 4][k])/d11;
            q[ 6*length+k] = (U*p[ 9][k] + v*p[ 8][k])/d12;
            q[ 9*length+k] = (u*p[19][k] + V*p[18][k])/d21;
            q[10*length+k] = (U*p[13][k] + V*p[14][k])/d22;
        }

        memcpy(q+ 0*length, p[ 0], length*sizeof(float));
        memcpy(q+ 1*length, p[ 1], length*sizeof(float));
        memcpy(q+ 2*length, p[ 7], length*sizeof(float));
        memcpy(q+ 3*length, p[ 5], length*sizeof(float));
        memcpy(q+ 4*length, p[ 2], length*sizeof(float));
        memcpy(q+ 7*length, p[ 6], length*sizeof(float));
        memcpy(q+ 8*length, p[16], length*sizeof(float));
        memcpy(q+11*length, p[12], length*sizeof(float));
        memcpy(q+12*length, p[15], length*sizeof(float));
        memcpy(q+13*length, p[17], length*sizeof(float));
        memcpy(q+14*length, p[11], length*sizeof(float));
        memcpy(q+15*length, p[10], length*sizeof(float));

        if (badPatch) {
            memset(q, 0, 16*3*sizeof(float));
        }

        for (int x = 0; x < 4; ++x) {
            for (int y = 0; y < 4; ++y) {
                int index = y*4 + x;
                min[0] = std::min(min[0], q[index*3+0]);
                min[1] = std::min(min[1], q[index*3+1]);
                min[2] = std::min(min[2], q[index*3+2]);
                max[0] = std::max(max[0], q[index*3+0]);
                max[1] = std::max(max[1], q[index*3+1]);
                max[2] = std::max(max[2], q[index*3+2]);
                bezierVertices.push_back(q[index*3+0]);
                bezierVertices.push_back(q[index*3+1]);
                bezierVertices.push_back(q[index*3+2]);
            }
        }
        bezierBounds.push_back(min[0]);
        bezierBounds.push_back(min[1]);
        bezierBounds.push_back(min[2]);
        bezierBounds.push_back(max[0]);
        bezierBounds.push_back(max[1]);
        bezierBounds.push_back(max[2]);

#if 0
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 0]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 1]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 2]);
        cpIndices.push_back(patchTables->GetPatchTable()[parray.GetVertIndex() + patchIndex*4 + 3]);
#else
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
        cpIndices.push_back(-1);
#endif
    }

    return parray.GetNumPatches();
}
