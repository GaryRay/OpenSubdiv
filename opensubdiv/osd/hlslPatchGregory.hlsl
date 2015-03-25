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

#if defined OSD_FRACTIONAL_ODD_SPACING
    #define HS_PARTITION "fractional_odd"
#elif defined OSD_FRACTIONAL_EVEN_SPACING
    #define HS_PARTITION "fractional_even"
#else
    #define HS_PARTITION "integer"
#endif

//----------------------------------------------------------
// Patches.Coefficients
//----------------------------------------------------------

#if OSD_MAX_VALENCE<=10
static float ef[7] = {
    0.813008, 0.500000, 0.363636, 0.287505,
    0.238692, 0.204549, 0.179211
};
#else
static float ef[27] = {
    0.812816, 0.500000, 0.363644, 0.287514,
    0.238688, 0.204544, 0.179229, 0.159657,
    0.144042, 0.131276, 0.120632, 0.111614,
    0.103872, 0.09715, 0.0912559, 0.0860444,
    0.0814022, 0.0772401, 0.0734867, 0.0700842,
    0.0669851, 0.0641504, 0.0615475, 0.0591488,
    0.0569311, 0.0548745, 0.0529621
};
#endif

float cosfn(uint n, uint j) {
    return cos((2.0f * M_PI * j)/float(n));
}

float sinfn(uint n, uint j) {
    return sin((2.0f * M_PI * j)/float(n));    
}


//----------------------------------------------------------
// Patches.TessVertexGregory
//----------------------------------------------------------

Buffer<float> VertexBuffer : register( t0 );
Buffer<int> OsdValenceBuffer : register( t1 );

void vs_main_patches( in InputVertex input,
                      uint vID : SV_VertexID,
                      out GregHullVertex output )
{
    output.hullPosition = mul(OsdModelViewMatrix(), input.position).xyz;
    OSD_PATCH_CULL_COMPUTE_CLIPFLAGS(input.position);

    int ivalence = OsdValenceBuffer[int(vID * (2 * OSD_MAX_VALENCE + 1))];
    output.valence = ivalence;
    uint valence = uint(abs(ivalence));

    float3 f[OSD_MAX_VALENCE];
    float3 pos = input.position.xyz;
    float3 opos = float3(0,0,0);

#ifdef OSD_PATCH_GREGORY_BOUNDARY
    output.org = input.position.xyz;
    int boundaryEdgeNeighbors[2];
    uint currNeighbor = 0;
    uint ibefore = 0;
    uint zerothNeighbor = 0;
#endif

    for (uint i=0; i<valence; ++i) {
        uint im=(i+valence-1)%valence;
        uint ip=(i+1)%valence;

        uint idx_neighbor = uint(OsdValenceBuffer[int(vID * (2*OSD_MAX_VALENCE+1) + 2*i + 0 + 1)]);

#ifdef OSD_PATCH_GREGORY_BOUNDARY
        bool isBoundaryNeighbor = false;
        int valenceNeighbor = OsdValenceBuffer[int(idx_neighbor * (2*OSD_MAX_VALENCE+1))];

        if (valenceNeighbor < 0) {
            isBoundaryNeighbor = true;
            if (currNeighbor<2) {
                boundaryEdgeNeighbors[currNeighbor] = int(idx_neighbor);
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
#endif

        float3 neighbor =
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor+2)]);

        uint idx_diagonal = uint(OsdValenceBuffer[int(vID * (2*OSD_MAX_VALENCE+1) + 2*i + 1 + 1)]);

        float3 diagonal =
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal+2)]);

        uint idx_neighbor_p = uint(OsdValenceBuffer[int(vID * (2*OSD_MAX_VALENCE+1) + 2*ip + 0 + 1)]);

        float3 neighbor_p =
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor_p)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor_p+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor_p+2)]);

        uint idx_neighbor_m = uint(OsdValenceBuffer[int(vID * (2*OSD_MAX_VALENCE+1) + 2*im + 0 + 1)]);

        float3 neighbor_m =
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor_m)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor_m+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor_m+2)]);

        uint idx_diagonal_m = uint(OsdValenceBuffer[int(vID * (2*OSD_MAX_VALENCE+1) + 2*im + 1 + 1)]);

        float3 diagonal_m =
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal_m)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal_m+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal_m+2)]);

        f[i] = (pos * float(valence) + (neighbor_p + neighbor)*2.0f + diagonal) / (float(valence)+5.0f);

        opos += f[i];
        output.r[i] = (neighbor_p-neighbor_m)/3.0f + (diagonal - diagonal_m)/6.0f;
    }

    opos /= valence;
    output.position = float4(opos, 1.0f).xyz;

    float3 e;
    output.e0 = float3(0,0,0);
    output.e1 = float3(0,0,0);

    for(uint i=0; i<valence; ++i) {
        uint im = (i + valence -1) % valence;
        e = 0.5f * (f[i] + f[im]);
        output.e0 += cosfn(valence, i)*e;
        output.e1 += sinfn(valence, i)*e;
    }
    output.e0 *= ef[valence - 3];
    output.e1 *= ef[valence - 3];

#ifdef OSD_PATCH_GREGORY_BOUNDARY
    output.zerothNeighbor = zerothNeighbor;
    if (currNeighbor == 1) {
        boundaryEdgeNeighbors[1] = boundaryEdgeNeighbors[0];
    }

    if (ivalence < 0) {
        if (valence > 2) {
            output.position = (
                float3(VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0])],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0]+1)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0]+2)]) +
                float3(VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1])],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1]+1)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1]+2)]) +
                4.0f * pos)/6.0f;
        } else {
            output.position = pos;
        }

        output.e0 = (
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0])],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0]+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0]+2)]) -
            float3(VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1])],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1]+1)],
                   VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1]+2)])
            )/6.0;

        float k = float(float(valence) - 1.0f);    //k is the number of faces
        float c = cos(M_PI/k);
        float s = sin(M_PI/k);
        float gamma = -(4.0f*s)/(3.0f*k+c);
        float alpha_0k = -((1.0f+2.0f*c)*sqrt(1.0f+c))/((3.0f*k+c)*sqrt(1.0f-c));
        float beta_0 = s/(3.0f*k + c);


        int idx_diagonal = OsdValenceBuffer[int((vID) * (2*OSD_MAX_VALENCE+1) + 2*zerothNeighbor + 1 + 1)];
        idx_diagonal = abs(idx_diagonal);
        float3 diagonal =
                float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal+1)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal+2)]);

        output.e1 = gamma * pos +
            alpha_0k * float3(VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0])],
                              VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0]+1)],
                              VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[0]+2)]) +
            alpha_0k * float3(VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1])],
                              VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1]+1)],
                              VertexBuffer[int(OSD_NUM_ELEMENTS*boundaryEdgeNeighbors[1]+2)]) +
            beta_0 * diagonal;

        for (uint x=1; x<valence - 1; ++x) {
            uint curri = ((x + zerothNeighbor)%valence);
            float alpha = (4.0f*sin((M_PI * float(x))/k))/(3.0f*k+c);
            float beta = (sin((M_PI * float(x))/k) + sin((M_PI * float(x+1))/k))/(3.0f*k+c);

            int idx_neighbor = OsdValenceBuffer[int((vID) * (2*OSD_MAX_VALENCE+1) + 2*curri + 0 + 1)];
            idx_neighbor = abs(idx_neighbor);

            float3 neighbor =
                float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor+1)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*idx_neighbor+2)]);

            idx_diagonal = OsdValenceBuffer[int((vID) * (2*OSD_MAX_VALENCE+1) + 2*curri + 1 + 1)];

            diagonal =
                float3(VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal+1)],
                       VertexBuffer[int(OSD_NUM_ELEMENTS*idx_diagonal+2)]);

            output.e1 += alpha * neighbor + beta * diagonal;
        }

        output.e1 /= 3.0f;
    }
#endif
}

//----------------------------------------------------------
// Patches.HullGregory
//----------------------------------------------------------

Buffer<int> OsdQuadOffsetBuffer : register( t2 );

HS_CONSTANT_FUNC_OUT HSConstFunc(
    InputPatch<GregHullVertex, 4> patch,
    uint primitiveID : SV_PrimitiveID)
{
    HS_CONSTANT_FUNC_OUT output;
    int patchLevel = GetPatchLevel(primitiveID);

    OSD_PATCH_CULL(4);

#ifdef OSD_ENABLE_SCREENSPACE_TESSELLATION
    output.tessLevelOuter[0] =
        TessAdaptive(patch[0].hullPosition.xyz, patch[1].hullPosition.xyz);
    output.tessLevelOuter[1] =
        TessAdaptive(patch[0].hullPosition.xyz, patch[3].hullPosition.xyz);
    output.tessLevelOuter[2] =
        TessAdaptive(patch[2].hullPosition.xyz, patch[3].hullPosition.xyz);
    output.tessLevelOuter[3] =
        TessAdaptive(patch[1].hullPosition.xyz, patch[2].hullPosition.xyz);
    output.tessLevelInner[0] =
        max(output.tessLevelOuter[1], output.tessLevelOuter[3]);
    output.tessLevelInner[1] =
        max(output.tessLevelOuter[0], output.tessLevelOuter[2]);
#else
    output.tessLevelInner[0] = GetTessLevel(patchLevel);
    output.tessLevelInner[1] = GetTessLevel(patchLevel);
    output.tessLevelOuter[0] = GetTessLevel(patchLevel);
    output.tessLevelOuter[1] = GetTessLevel(patchLevel);
    output.tessLevelOuter[2] = GetTessLevel(patchLevel);
    output.tessLevelOuter[3] = GetTessLevel(patchLevel);
#endif
    return output;
}

[domain("quad")]
[partitioning(HS_PARTITION)]
[outputtopology("triangle_ccw")]
[outputcontrolpoints(4)]
[patchconstantfunc("HSConstFunc")]
GregDomainVertex hs_main_patches(
    in InputPatch<GregHullVertex, 4> patch,
    uint primitiveID : SV_PrimitiveID,
    in uint ID : SV_OutputControlPointID )
{
    uint i = ID;
    uint ip = (i+1)%4;
    uint im = (i+3)%4;
    uint valence = abs(patch[i].valence);
    uint n = valence;
    int base = OsdGregoryQuadOffsetBase();

    GregDomainVertex output;
    output.position = patch[ID].position;

    uint start = uint(OsdQuadOffsetBuffer[int(4*primitiveID+base + i)]) & 0x00ffu;
    uint prev = uint(OsdQuadOffsetBuffer[int(4*primitiveID+base + i)]) & 0xff00u;
    prev = uint(prev/256);

    uint start_m = uint(OsdQuadOffsetBuffer[int(4*primitiveID+base + im)]) & 0x00ffu;
    uint prev_p = uint(OsdQuadOffsetBuffer[int(4*primitiveID+base + ip)]) & 0xff00u;
    prev_p = uint(prev_p/256);

    uint np = abs(patch[ip].valence);
    uint nm = abs(patch[im].valence);

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

#ifdef OSD_PATCH_GREGORY_BOUNDARY
    float3 Ep = float3(0.0f,0.0f,0.0f);
    float3 Em = float3(0.0f,0.0f,0.0f);
    float3 Fp = float3(0.0f,0.0f,0.0f);
    float3 Fm = float3(0.0f,0.0f,0.0f);

    float3 Em_ip;
    if (patch[ip].valence < -2) {
        uint j = (np + prev_p - patch[ip].zerothNeighbor) % np;
        Em_ip = patch[ip].position + cos((M_PI*j)/float(np-1))*patch[ip].e0 + sin((M_PI*j)/float(np-1))*patch[ip].e1;
    } else {
        Em_ip = patch[ip].position + patch[ip].e0*cosfn(np, prev_p) + patch[ip].e1*sinfn(np, prev_p);
    }

    float3 Ep_im;
    if (patch[im].valence < -2) {
        uint j = (nm + start_m - patch[im].zerothNeighbor) % nm;
        Ep_im = patch[im].position + cos((M_PI*j)/float(nm-1))*patch[im].e0 + sin((M_PI*j)/float(nm-1))*patch[im].e1;
    } else {
        Ep_im = patch[im].position + patch[im].e0*cosfn(nm, start_m) + patch[im].e1*sinfn(nm, start_m);
    }

    if (patch[i].valence < 0) {
        n = (n-1)*2;
    }
    if (patch[im].valence < 0) {
        nm = (nm-1)*2;
    }
    if (patch[ip].valence < 0) {
        np = (np-1)*2;
    }

    if (patch[i].valence > 2) {
        Ep = patch[i].position + (patch[i].e0*cosfn(n, start) + patch[i].e1*sinfn(n, start));
        Em = patch[i].position + (patch[i].e0*cosfn(n, prev) +  patch[i].e1*sinfn(n, prev));

        float s1=3-2*cosfn(n,1)-cosfn(np,1);
        float s2=2*cosfn(n,1);

        Fp = (cosfn(np,1)*patch[i].position + s1*Ep + s2*Em_ip + patch[i].r[start])/3.0f;
        s1 = 3.0f-2.0f*cos(2.0f*M_PI/float(n))-cos(2.0f*M_PI/float(nm));
        Fm = (cosfn(nm,1)*patch[i].position + s1*Em + s2*Ep_im - patch[i].r[prev])/3.0f;

    } else if (patch[i].valence < -2) {
        uint j = (valence + start - patch[i].zerothNeighbor) % valence;

        Ep = patch[i].position + cos((M_PI*j)/float(valence-1))*patch[i].e0 + sin((M_PI*j)/float(valence-1))*patch[i].e1;
        j = (valence + prev - patch[i].zerothNeighbor) % valence;
        Em = patch[i].position + cos((M_PI*j)/float(valence-1))*patch[i].e0 + sin((M_PI*j)/float(valence-1))*patch[i].e1;

        float3 Rp = ((-2.0f * patch[i].org - 1.0f * patch[im].org) + (2.0f * patch[ip].org + 1.0f * patch[(i+2)%4].org))/3.0f;
        float3 Rm = ((-2.0f * patch[i].org - 1.0f * patch[ip].org) + (2.0f * patch[im].org + 1.0f * patch[(i+2)%4].org))/3.0f;

        float s1 = 3-2*cosfn(n,1)-cosfn(np,1);
        float s2 = 2*cosfn(n,1);

        Fp = (cosfn(np,1)*patch[i].position + s1*Ep + s2*Em_ip + patch[i].r[start])/3.0f;
        s1 = 3.0f-2.0f*cos(2.0f*M_PI/float(n))-cos(2.0f*M_PI/float(nm));
        Fm = (cosfn(nm,1)*patch[i].position + s1*Em + s2*Ep_im - patch[i].r[prev])/3.0f;

        if (patch[im].valence < 0) {
            s1=3-2*cosfn(n,1)-cosfn(np,1);
            Fp = Fm = (cosfn(np,1)*patch[i].position + s1*Ep + s2*Em_ip + patch[i].r[start])/3.0f;
        } else if (patch[ip].valence < 0) {
            s1 = 3.0f-2.0f*cos(2.0f*M_PI/n)-cos(2.0f*M_PI/nm);
            Fm = Fp = (cosfn(nm,1)*patch[i].position + s1*Em + s2*Ep_im - patch[i].r[prev])/3.0f;
        }

    } else if (patch[i].valence == -2) {
        Ep = (2.0f * patch[i].org + patch[ip].org)/3.0f;
        Em = (2.0f * patch[i].org + patch[im].org)/3.0f;
        Fp = Fm = (4.0f * patch[i].org + patch[(i+2)%n].org + 2.0f * patch[ip].org + 2.0f * patch[im].org)/9.0f;
    }

#else // not OSD_PATCH_GREGORY_BOUNDARY

    float3 Ep = patch[i].position + patch[i].e0 * cosfn(n, start) + patch[i].e1*sinfn(n, start);
    float3 Em = patch[i].position + patch[i].e0 * cosfn(n, prev ) + patch[i].e1*sinfn(n, prev );

    float3 Em_ip = patch[ip].position + patch[ip].e0*cosfn(np, prev_p) + patch[ip].e1*sinfn(np, prev_p);
    float3 Ep_im = patch[im].position + patch[im].e0*cosfn(nm, start_m) + patch[im].e1*sinfn(nm, start_m);

    float s1 = 3-2*cosfn(n,1)-cosfn(np,1);
    float s2 = 2*cosfn(n,1);

    float3 Fp = (cosfn(np,1)*patch[i].position + s1*Ep + s2*Em_ip + patch[i].r[start])/3.0f;
    s1 = 3.0f-2.0f*cos(2.0f*M_PI/float(n))-cos(2.0f*M_PI/float(nm));
    float3 Fm = (cosfn(nm,1)*patch[i].position + s1*Em +s2*Ep_im - patch[i].r[prev])/3.0f;

#endif

    output.Ep = Ep;
    output.Em = Em;
    output.Fp = Fp;
    output.Fm = Fm;

    int patchLevel = GetPatchLevel(primitiveID);
    output.patchCoord = float4(0, 0,
                               patchLevel+0.5f,
                               GetPrimitiveID(primitiveID)+0.5f);

    OSD_COMPUTE_PTEX_COORD_HULL_SHADER;

    return output;
}

//----------------------------------------------------------
// Patches.DomainGregory
//----------------------------------------------------------

[domain("quad")]
void ds_main_patches(
    in HS_CONSTANT_FUNC_OUT input,
    in OutputPatch<GregDomainVertex, 4> patch,
    in float2 uv : SV_DomainLocation,
    out OutputVertex output )
{
    float u = uv.x,
          v = uv.y;

    float3 p[20];

    p[0] = patch[0].position;
    p[1] = patch[0].Ep;
    p[2] = patch[0].Em;
    p[3] = patch[0].Fp;
    p[4] = patch[0].Fm;

    p[5] = patch[1].position;
    p[6] = patch[1].Ep;
    p[7] = patch[1].Em;
    p[8] = patch[1].Fp;
    p[9] = patch[1].Fm;

    p[10] = patch[2].position;
    p[11] = patch[2].Ep;
    p[12] = patch[2].Em;
    p[13] = patch[2].Fp;
    p[14] = patch[2].Fm;

    p[15] = patch[3].position;
    p[16] = patch[3].Ep;
    p[17] = patch[3].Em;
    p[18] = patch[3].Fp;
    p[19] = patch[3].Fm;

    float3 q[16];

    float U = 1-u, V=1-v;

    float d11 = u+v; if(u+v==0.0f) d11 = 1.0f;
    float d12 = U+v; if(U+v==0.0f) d12 = 1.0f;
    float d21 = u+V; if(u+V==0.0f) d21 = 1.0f;
    float d22 = U+V; if(U+V==0.0f) d22 = 1.0f;

    q[ 5] = (u*p[3] + v*p[4])/d11;
    q[ 6] = (U*p[9] + v*p[8])/d12;
    q[ 9] = (u*p[19] + V*p[18])/d21;
    q[10] = (U*p[13] + V*p[14])/d22;

    q[ 0] = p[0];
    q[ 1] = p[1];
    q[ 2] = p[7];
    q[ 3] = p[5];
    q[ 4] = p[2];
    q[ 7] = p[6];
    q[ 8] = p[16];
    q[11] = p[12];
    q[12] = p[15];
    q[13] = p[17];
    q[14] = p[11];
    q[15] = p[10];

    float3 WorldPos  = float3(0, 0, 0);
    float3 Tangent   = float3(0, 0, 0);
    float3 BiTangent = float3(0, 0, 0);

#line 519

#ifdef OSD_COMPUTE_NORMAL_DERIVATIVES
    float B[4], D[4], C[4];

    float3 BUCP[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)},
           DUCP[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)},
           CUCP[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    float3 dUU = float3(0, 0, 0);
    float3 dVV = float3(0, 0, 0);
    float3 dUV = float3(0, 0, 0);

    Univar4x4(u, B, D, C);

    for (int i=0; i<4; ++i) {
        for (uint j=0; j<4; ++j) {
            // reverse face front
            float3 A = q[i + 4*j];

            BUCP[i] += A * B[j];
            DUCP[i] += A * D[j];
            CUCP[i] += A * C[j];
        }
    }

    Univar4x4(v, B, D, C);

    for (int i=0; i<4; ++i) {
        WorldPos  += B[i] * BUCP[i];
        Tangent   += B[i] * DUCP[i];
        BiTangent += D[i] * BUCP[i];
        dUU += B[i] * CUCP[i];
        dVV += C[i] * BUCP[i];
        dUV += D[i] * DUCP[i];
    }

    int level = int(patch[0].ptexInfo.z);
    BiTangent *= 3 * level;
    Tangent *= 3 * level;
    dUU *= 6 * level;
    dVV *= 6 * level;
    dUV *= 9 * level;

    float3 n = cross(Tangent, BiTangent);
    float3 normal = normalize(n);

    float E = dot(Tangent, Tangent);
    float F = dot(Tangent, BiTangent);
    float G = dot(BiTangent, BiTangent);
    float e = dot(normal, dUU);
    float f = dot(normal, dUV);
    float g = dot(normal, dVV);

    float3 Nu = (f*F-e*G)/(E*G-F*F) * Tangent + (e*F-f*E)/(E*G-F*F) * BiTangent;
    float3 Nv = (g*F-f*G)/(E*G-F*F) * Tangent + (f*F-g*E)/(E*G-F*F) * BiTangent;

    Nu = Nu/length(n) - n * (dot(Nu,n)/pow(dot(n,n), 1.5));
    Nv = Nv/length(n) - n * (dot(Nv,n)/pow(dot(n,n), 1.5));

    BiTangent = mul(OsdModelViewMatrix(), float4(BiTangent, 0)).xyz;
    Tangent = mul(OsdModelViewMatrix(), float4(Tangent, 0)).xyz;

    normal = normalize(cross(BiTangent, Tangent));

    output.Nu = Nu;
    output.Nv = Nv;

#else
    float B[4], D[4];
    float3 BUCP[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)},
           DUCP[4] = {float3(0,0,0), float3(0,0,0), float3(0,0,0), float3(0,0,0)};

    Univar4x4(uv.x, B, D);

    for (int i=0; i<4; ++i) {
        for (uint j=0; j<4; ++j) {
            // reverse face front
            float3 A = q[i + 4*j];

            BUCP[i] += A * B[j];
            DUCP[i] += A * D[j];
        }
    }

    Univar4x4(uv.y, B, D);

    for (uint i=0; i<4; ++i) {
        WorldPos  += B[i] * BUCP[i];
        Tangent   += B[i] * DUCP[i];
        BiTangent += D[i] * BUCP[i];
    }
    int level = int(patch[0].ptexInfo.z);
    BiTangent *= 3 * level;
    Tangent *= 3 * level;

    BiTangent = mul(OsdModelViewMatrix(), float4(BiTangent, 0)).xyz;
    Tangent = mul(OsdModelViewMatrix(), float4(Tangent, 0)).xyz;

    float3 normal = normalize(cross(BiTangent, Tangent));

#endif

    output.position = mul(OsdModelViewMatrix(), float4(WorldPos, 1.0f));
    output.normal = normal;
    output.tangent = BiTangent;
    output.bitangent = Tangent;

    output.patchCoord = patch[0].patchCoord;
    output.patchCoord.xy = float2(v, u);
	
	output.edgeDistance = 0;

    OSD_COMPUTE_PTEX_COORD_DOMAIN_SHADER;

    OSD_DISPLACEMENT_CALLBACK;

    output.positionOut = mul(OsdProjectionMatrix(),
                             float4(output.position.xyz, 1.0f));
}
