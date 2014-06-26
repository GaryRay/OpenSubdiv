struct BVHNode
{
    float bmin[3];
    float bmax[3];
    int flag; // 1 = leaf node, 0 = branch node
    int axis;

    // leaf
    //   data[0] = npoints
    //   data[1] = index
    //
    // branch
    //   data[0] = child[0]
    //   data[1] = child[1]
    unsigned int data[2];
};

struct Intersection {
    float t, u, v;
    unsigned int faceID;
    float4 normal;
    int clipLevel;
};

struct Ray {
    float4 org;
    float4 dir;
    float4 invDir;
    int dirSign[4]; // dirSign[3] is used for the destination index of the resulting image.
};

struct BezierPatch {
    float4 cp[16];
};

struct BezierPatchIntersection {
    struct BezierPatch patch;
    float4 min;
    float4 max;
    float eps;
    float uRange[2];
    float vRange[2];
};

struct RangeAABB {
    float tmin, tmax;
};

struct UVT {
    float u, v, t;
    int clipLevel;
};

struct Matrix4 {
    float4 v[4];
};

// --------------------------------------------------------------------
// prototypes
bool IntersectRayAABB(float *tminOut, float *tmaxOut,
                      float maxT, float *bmin, float *bmax,
                      float4 rayOrg, float4 rayInvDir, int *rayDirSign);
bool TestLeafNode(struct Intersection *isect, // [inout]
                  const struct BVHNode *node,
                  __global const unsigned int *indices,
                  __global const float *bezierVerts,
                  const struct Ray *ray);
bool PatchIsect(struct Intersection *isect,
                __global const float *bezierVerts,
                const struct Ray *ray);
bool TriangleIsect(float *tInOut, float *uOut, float *vOut,
                   float4 p0, float4 p1, float4 p2,
                   float4 rayOrg, float4 rayDir);

void BpConstruct(struct BezierPatch *patch, __global const float *verts);
void BpGetMinMax(const struct BezierPatch *patch, float4 *min, float4 *max, float eps);
void BpTransform(struct BezierPatch *patch, const struct Matrix4 *mat);
void BpRotate(struct BezierPatch *patch, float4 dx);
void BpRotateU(struct BezierPatch *patch);
void BpRotateU(struct BezierPatch *patch);
float4 BpEvaluateDu(struct BezierPatch *patch, float u, float v);
float4 BpEvaluateDv(struct BezierPatch *patch, float u, float v);

void BpiConstruct(struct BezierPatchIntersection *isect,
                  const struct BezierPatch *patch);
bool BpiTest(struct BezierPatchIntersection *isect,
             struct Intersection *info,
             const struct Ray *r,
             float tmin, float tmax);
void BezierSplit(float4 a[], float4 b[], const float4 c[], float t, int stride);

bool IntersectAABB(struct RangeAABB *rng, float4 min, float4 max,
                   const struct Ray *r, float tmin, float tmax);
bool BpiTestInternal(struct BezierPatchIntersection *isect,
                     struct Intersection *info, const struct Ray *r,
                     float tmin, float tmax);
bool BpiTestBezierPatch(struct BezierPatchIntersection *this,
                        struct UVT *uvt, struct BezierPatch *patch,
                        float tmin, float tmax, float eps);
bool BpiTestBezierClipL(struct BezierPatchIntersection *this,
                        struct UVT* info, struct BezierPatch *patch,
                        float u0, float u1, float v0, float v1,
                        float zmin, float zmax, int level);

void MatrixMultiply(struct Matrix4 *dst, const struct Matrix4 *a, const struct Matrix4 *b);
float4 MatrixApply(const struct Matrix4 *m, float4 v);
void GetZAlign(struct Matrix4 *mat, const struct Ray *r);
// --------------------------------------------------------------------

void MatrixMultiply(struct Matrix4 *dst, const struct Matrix4 *lhs, const struct Matrix4 *rhs)
{
    dst->v[0] = (float4)(
        lhs->v[0].x*rhs->v[0].x + lhs->v[0].y*rhs->v[1].x + lhs->v[0].z*rhs->v[2].x + lhs->v[0].w*rhs->v[3].x,
        lhs->v[0].x*rhs->v[0].y + lhs->v[0].y*rhs->v[1].y + lhs->v[0].z*rhs->v[2].y + lhs->v[0].w*rhs->v[3].y,
        lhs->v[0].x*rhs->v[0].z + lhs->v[0].y*rhs->v[1].z + lhs->v[0].z*rhs->v[2].z + lhs->v[0].w*rhs->v[3].z,
        lhs->v[0].x*rhs->v[0].w + lhs->v[0].y*rhs->v[1].w + lhs->v[0].z*rhs->v[2].w + lhs->v[0].w*rhs->v[3].w );

    dst->v[1] = (float4)(
        lhs->v[1].x*rhs->v[0].x + lhs->v[1].y*rhs->v[1].x + lhs->v[1].z*rhs->v[2].x + lhs->v[1].w*rhs->v[3].x,
        lhs->v[1].x*rhs->v[0].y + lhs->v[1].y*rhs->v[1].y + lhs->v[1].z*rhs->v[2].y + lhs->v[1].w*rhs->v[3].y,
        lhs->v[1].x*rhs->v[0].z + lhs->v[1].y*rhs->v[1].z + lhs->v[1].z*rhs->v[2].z + lhs->v[1].w*rhs->v[3].z,
        lhs->v[1].x*rhs->v[0].w + lhs->v[1].y*rhs->v[1].w + lhs->v[1].z*rhs->v[2].w + lhs->v[1].w*rhs->v[3].w );

    dst->v[2] = (float4)(
        lhs->v[2].x*rhs->v[0].x + lhs->v[2].y*rhs->v[1].x + lhs->v[2].z*rhs->v[2].x + lhs->v[2].w*rhs->v[3].x,
        lhs->v[2].x*rhs->v[0].y + lhs->v[2].y*rhs->v[1].y + lhs->v[2].z*rhs->v[2].y + lhs->v[2].w*rhs->v[3].y,
        lhs->v[2].x*rhs->v[0].z + lhs->v[2].y*rhs->v[1].z + lhs->v[2].z*rhs->v[2].z + lhs->v[2].w*rhs->v[3].z,
        lhs->v[2].x*rhs->v[0].w + lhs->v[2].y*rhs->v[1].w + lhs->v[2].z*rhs->v[2].w + lhs->v[2].w*rhs->v[3].w );

    dst->v[3] = (float4)(
        lhs->v[3].x*rhs->v[0].x + lhs->v[3].y*rhs->v[1].x + lhs->v[3].z*rhs->v[2].x + lhs->v[3].w*rhs->v[3].x,
        lhs->v[3].x*rhs->v[0].y + lhs->v[3].y*rhs->v[1].y + lhs->v[3].z*rhs->v[2].y + lhs->v[3].w*rhs->v[3].y,
        lhs->v[3].x*rhs->v[0].z + lhs->v[3].y*rhs->v[1].z + lhs->v[3].z*rhs->v[2].z + lhs->v[3].w*rhs->v[3].z,
        lhs->v[3].x*rhs->v[0].w + lhs->v[3].y*rhs->v[1].w + lhs->v[3].z*rhs->v[2].w + lhs->v[3].w*rhs->v[3].w );
}

float4 MatrixApply(const struct Matrix4 *m, float4 v)
{
    float4 r;
    r.x = m->v[0].x * v.x + m->v[0].y * v.y + m->v[0].z * v.z + m->v[0].w;
    r.y = m->v[1].x * v.x + m->v[1].y * v.y + m->v[1].z * v.z + m->v[1].w;
    r.z = m->v[2].x * v.x + m->v[2].y * v.y + m->v[2].z * v.z + m->v[2].w;
    r.w = m->v[3].x * v.x + m->v[3].y * v.y + m->v[3].z * v.z + m->v[3].w;
    float ir = 1.0f/r.w;
    r.x *= ir;
    r.y *= ir;
    r.z *= ir;
    r.w = 0;
    return r;
}

void GetZAlign(struct Matrix4 *mat, const struct Ray *r)
{
    float4 org = r->org;
    float4 dir = r->dir;
    float z[4] = {dir.x, dir.y, dir.z, 0};

    int plane = 0;
    if (fabs(z[1]) < fabs(z[plane])) plane = 1;
    if (fabs(z[2]) < fabs(z[plane])) plane = 2;

    float4 x = (float4)(0, 0, 0, 0);
    if (plane == 0) x.x = 1;
    if (plane == 1) x.y = 1;
    if (plane == 2) x.z = 1;
    float4 y = normalize(cross(dir,x));
    x = cross(y,dir);

    x.w = y.w = dir.w = 0;

    struct Matrix4 rot;
    rot.v[0] = x;
    rot.v[1] = y;
    rot.v[2] = dir;
    rot.v[3] = (float4)(0, 0, 0, 1);

    struct Matrix4 trs;
    trs.v[0] = (float4)(1, 0, 0, -org.x);
    trs.v[1] = (float4)(0, 1, 0, -org.y);
    trs.v[2] = (float4)(0, 0, 1, -org.z);
    trs.v[3] = (float4)(0, 0, 0, 1);

    MatrixMultiply(mat, &rot, &trs);
}

void BpConstruct(struct BezierPatch *patch, __global const float *verts)
{
#pragma unroll
    for (int i = 0; i < 16; ++i) {
        patch->cp[i] = (float4)(verts[i*3+0], verts[i*3+1], verts[i*3+2], 0);
    }
}
void BpGetMinMax(const struct BezierPatch *patch, float4 *min, float4 *max, float eps)
{
    *min = patch->cp[0];
    *max = patch->cp[0];
    for (int i = 1; i < 16; ++i) {
        *min = fmin(*min, patch->cp[i]);
        *max = fmax(*max, patch->cp[i]);
    }
    *min -= (float4)(eps, eps, eps, 0);
    *max += (float4)(eps, eps, eps, 0);
}

void BpTransform(struct BezierPatch *patch, const struct Matrix4 *mat)
{
    for (int i = 0; i < 16; ++i) {
        patch->cp[i] = MatrixApply(mat, patch->cp[i]);
    }
}

void BpiConstruct(struct BezierPatchIntersection *isect,
                  const struct BezierPatch *patch)
{
    BpGetMinMax(patch, &isect->min, &isect->max, 0.01f);
    isect->patch = *patch;
    isect->eps = 1.0e-4f;
}

bool IntersectAABB(struct RangeAABB *rng, float4 min, float4 max,
                   const struct Ray *r, float tmin, float tmax)
{
    float4 box[2];
    box[0] = min;
    box[1] = max;
    int4 sign = (int4)(r->dirSign[0], r->dirSign[1], r->dirSign[2], 0);

    float4 org = r->org;
    float4 idir = r->invDir;

    tmin = fmax(tmin, (box[  sign.x].x-org.x)*idir.x);
    tmin = fmax(tmin, (box[  sign.y].y-org.y)*idir.y);
    tmin = fmax(tmin, (box[  sign.z].z-org.z)*idir.z);
    tmax = fmin(tmax, (box[1-sign.x].x-org.x)*idir.x);
    tmax = fmin(tmax, (box[1-sign.y].y-org.y)*idir.y);
    tmax = fmin(tmax, (box[1-sign.z].z-org.z)*idir.z);

    rng->tmin = tmin;
    rng->tmax = tmax;

    return rng->tmin <= rng->tmax;
}

bool BpiTest(struct BezierPatchIntersection *isect,
             struct Intersection *info,
             const struct Ray *r,
             float tmin, float tmax)
{
    struct RangeAABB rng;
    if (IntersectAABB(&rng, isect->min, isect->max, r, tmin, tmax)) {
        tmin = fmax(tmin, rng.tmin);
        tmax = fmin(tmax, rng.tmax);

        return BpiTestInternal(isect, info, r, tmin, tmax);
    }
    return false;
}

bool BpiTestInternal(struct BezierPatchIntersection *isect,
                     struct Intersection *info, const struct Ray *r,
                     float tmin, float tmax)
{
    struct Matrix4 mat;
    GetZAlign(&mat, r);

    struct UVT uvt;
    uvt.t = tmax;
    struct BezierPatch patch = isect->patch;
    BpTransform(&patch, &mat);
    if (BpiTestBezierPatch(isect, &uvt, &patch, tmin, tmax, isect->eps)) {
        float t = uvt.t;
        float u = uvt.u;
        float v = uvt.v;

        info->t = t;
        info->u = u;
        info->v = v;
        info->clipLevel = uvt.clipLevel;
        {
            float4 du = BpEvaluateDu(&isect->patch, u, v);
            float4 dv = BpEvaluateDv(&isect->patch, u, v);
            info->normal = normalize(cross(du, dv));
        }
        return true;
    }
    return false;
}

float4 evaluateD(float t, float4 *cp)
{
    float t2 = t*t;
    return cp[0] * (3*t2*-1 + 2*t* 3 + -3)
         + cp[1] * (3*t2* 3 + 2*t*-6 +  3)
         + cp[2] * (3*t2*-3 + 2*t* 3)
         + cp[3] * (3*t2* 1);
}
float4 evaluate(float t, float4 *cp)
{
    return cp[0] * (1-t) * (1-t) * (1-t)
        + cp[1] * 3 * t * (1-t) * (1-t)
        + cp[2] * 3 * t * t * (1-t)
        + cp[3] * t * t * t;
}

float4 BpEvaluateDu(struct BezierPatch *patch, float u, float v)
{
    float4 b[4];
    for (int i = 0; i < 4; ++i) {
        b[i] = evaluateD(u, &patch->cp[i*4]);
    }
    return evaluate(v, b);
}

float4 BpEvaluateDv(struct BezierPatch *patch, float u, float v)
{
    float4 b[4];
    for (int i = 0; i < 4; ++i) {
        b[i] = evaluate(u, &patch->cp[i*4]);
    }
    return evaluateD(v, b);
}

void BpRotate(struct BezierPatch *patch, float4 dx)
{
    // normalize2
    float inv_len = 1.0f/sqrt(dx.x*dx.x + dx.y*dx.y);
    dx.x *= inv_len;
    dx.y *= inv_len;
    dx.z = 0;
    float4 dy = (float4)(-dx.y, dx.x, 0, 0);

    struct Matrix4 mat;
    mat.v[0] = dx;
    mat.v[1] = dy;
    mat.v[2] = (float4)(0, 0, 1, 0);
    mat.v[3] = (float4)(0, 0, 0, 1);

    BpTransform(patch, &mat);
}

void BpRotateU(struct BezierPatch *patch)
{
    float4 dx = patch->cp[12] - patch->cp[0] + patch->cp[15] - patch->cp[3];
    BpRotate(patch, dx);
}
void BpRotateV(struct BezierPatch *patch)
{
    float4 dx = patch->cp[3] - patch->cp[0] + patch->cp[15] - patch->cp[12];
    BpRotate(patch, dx);
}

void BezierSplit(float4 a[], float4 b[], const float4 c[], float t, int stride)
{
    float S = 1-t;
    float4 p0 = c[0*stride];
    float4 p1 = c[1*stride];
    float4 p2 = c[2*stride];
    float4 p3 = c[3*stride];
    a[0*stride] = p0;
    a[1*stride] = p0*S + p1*t;
    a[2*stride] = p0*S*S + p1*2*S*t + p2*t*t;
    a[3*stride] = p0*S*S*S + p1*3*S*S*t + p2*3*S*t*t + p3*t*t*t;

    b[0*stride] = p0*S*S*S + p1*3*S*S*t + p2*3*S*t*t + p3*t*t*t;
    b[1*stride] = p3*t*t + p2*2*t*S + p1*S*S;
    b[2*stride] = p3*t + p2*S;
    b[3*stride] = p3;
}

void BpSplitU(const struct BezierPatch *patch, struct BezierPatch patches[2], float u)
{
    for (int i = 0; i < 4; ++i) {
        BezierSplit(&patches[0].cp[i*4+0],
                    &patches[1].cp[i*4+0],
                    &patch->cp[i*4+0], u, 1);
    }
}
void BpSplitV(const struct BezierPatch *patch, struct BezierPatch patches[2], float v)
{
    for (int i = 0; i < 4; ++i) {
        BezierSplit(&patches[0].cp[i],
                    &patches[1].cp[i],
                    &patch->cp[i], v, 4);
    }
}

bool BpiTestBezierClipL(struct BezierPatchIntersection *this,
                        struct UVT* info, struct BezierPatch *patch,
                        float u0, float u1, float v0, float v1,
                        float zmin, float zmax, int level)
{
    // TODO.
    // DIRECT_BILINEAR
    float4 p0, p1, p2, p3;
    float4 rayOrg = (float4)(0, 0, 0, 0);
    float4 rayDir= (float4)(0, 0, 1, 0);
    p0 = patch->cp[0];
    p1 = patch->cp[3];
    p2 = patch->cp[12];
    p3 = patch->cp[15];
    bool bRet = false;
    float t = zmax, uu = 0, vv = 0;
    if (TriangleIsect(&t, &uu, &vv, p0, p2, p1, rayOrg, rayDir)) {
        float ww = 1 - (uu + vv);
        float u = ww*0 + uu*0 + vv*1;//00 - 01 - 10
        float v = ww*0 + uu*1 + vv*0;//00 - 01 - 10
        info->u = mix(u0,u1,u);
        info->v = mix(v0,v1,v);
        info->t = t;
        info->clipLevel = level;
        bRet = true;
    }
    if (TriangleIsect(&t, &uu, &vv, p1, p2, p3, rayOrg, rayDir)) {
        float ww = 1 - (uu + vv);
        float u = ww*1 + uu*0 + vv*1;//10 - 01 - 11
        float v = ww*0 + uu*1 + vv*1;//10 - 01 - 11
        info->u = mix(u0,u1,u);
        info->v = mix(v0,v1,v);
        info->t = t;
        info->clipLevel = level;
        bRet = true;
    }
    return bRet;
}

bool BpiTestBezierPatch(struct BezierPatchIntersection *this,
                        struct UVT *info, struct BezierPatch *patch,
                        float zmin, float zmax, float eps)
{
    float4 min, max;
    BpGetMinMax(patch, &min, &max, eps*1e-3f);

    if (0 < min.x || max.x < 0) return false;//x
    if (0 < min.y || max.y < 0) return false;//y
    if (max.z < zmin || zmax < min.z) return false;//z

    //return BpiTestBezierClipL(this, info, patch, 0, 1, 0, 1, zmin, zmax, 0);

    // non-recursive iteration
    bool bRet = false;
    struct Entry {
        struct BezierPatch patch;
        float u0, u1, v0, v1;
    };

#define MAX_STACK_DEPTH 100
    struct Entry patchStack[MAX_STACK_DEPTH];
    int stackIndex = 0;
    patchStack[stackIndex].patch = *patch;
    patchStack[stackIndex].u0 = 0;
    patchStack[stackIndex].u1 = 1;
    patchStack[stackIndex].v0 = 0;
    patchStack[stackIndex].v1 = 1;

    int maxLoop = 100;
    for (int i = 0; i < maxLoop && stackIndex >= 0; ++i) {

        // pop a patch
        float u0 = patchStack[stackIndex].u0;
        float u1 = patchStack[stackIndex].u1;
        float v0 = patchStack[stackIndex].v0;
        float v1 = patchStack[stackIndex].v1;
        struct BezierPatch currentPatch = patchStack[stackIndex].patch;
        --stackIndex;

        struct BezierPatch tpatch = currentPatch;
        float lenU = length(tpatch.cp[3] - tpatch.cp[0]);
        float lenV = length(tpatch.cp[12] - tpatch.cp[0]);
        bool clipU = lenU > lenV;

        // if it's small enough, test bilinear.
        if (lenU < 0.01f && lenV < 0.01f) {
            if(BpiTestBezierClipL(this, info, &currentPatch,
                                  u0, u1, v0, v1, zmin, zmax, i)) {
                // info is updated.
                zmax = info->t;
                bRet = true;
            }
            continue;
        }

        // rotate
        if (clipU) {
            BpRotateU(&tpatch);
        } else {
            BpRotateV(&tpatch);
        }
        float4 min, max;
        BpGetMinMax(&tpatch, &min, &max, eps*1e-3f);

        // out
        if (0 < min.x || max.x < 0 || 0 < min.y || max.y < 0 || max.z < zmin || zmax < min.z) {
            continue;
        }

        // push children
        struct Entry children[2];
        struct BezierPatch tmp[2];
        if (clipU) {
            BpSplitU(&currentPatch, tmp, 0.5f);
            float um = (u0+u1)*0.5f;
            float ut[4] = {u0,um,um,u1};

            for (int j = 0; j < 2; ++j) {
                ++stackIndex;
                patchStack[stackIndex].patch = tmp[j];
                patchStack[stackIndex].u0 = ut[2*j];
                patchStack[stackIndex].u1 = ut[2*j+1];
                patchStack[stackIndex].v0 = v0;
                patchStack[stackIndex].v1 = v1;
            }
        } else {
            BpSplitV(&currentPatch, tmp, 0.5f);
            float vm = (v0+v1)*0.5f;
            float vt[4] = {v0,vm,vm,v1};

            for (int j = 0; j < 2; ++j) {
                ++stackIndex;
                patchStack[stackIndex].patch = tmp[j];
                patchStack[stackIndex].u0 = u0;
                patchStack[stackIndex].u1 = u1;
                patchStack[stackIndex].v0 = vt[2*j];
                patchStack[stackIndex].v1 = vt[2*j+1];
            }
        }
        if (stackIndex >= MAX_STACK_DEPTH-1) break;
    }
    return bRet;
}

// --------------------------------------------------------------------

bool IntersectRayAABB(float *tminOut, float *tmaxOut,
                      float maxT, float *bmin, float *bmax,
                      float4 rayOrg, float4 rayInvDir, int *rayDirSign)
{
    float tmin, tmax;

    const float min_x = rayDirSign[0] * bmax[0] + (1-rayDirSign[0]) * bmin[0];
    const float min_y = rayDirSign[1] * bmax[1] + (1-rayDirSign[1]) * bmin[1];
    const float min_z = rayDirSign[2] * bmax[2] + (1-rayDirSign[2]) * bmin[2];
    const float max_x = rayDirSign[0] * bmin[0] + (1-rayDirSign[0]) * bmax[0];
    const float max_y = rayDirSign[1] * bmin[1] + (1-rayDirSign[1]) * bmax[1];
    const float max_z = rayDirSign[2] * bmin[2] + (1-rayDirSign[2]) * bmax[2];

    // X
    const float tmin_x = (min_x - rayOrg.x) * rayInvDir.x;
    const float tmax_x = (max_x - rayOrg.x) * rayInvDir.x;

    // Y
    const float tmin_y = (min_y - rayOrg.y) * rayInvDir.y;
    const float tmax_y = (max_y - rayOrg.y) * rayInvDir.y;

    tmin = (tmin_x > tmin_y) ? tmin_x : tmin_y;
    tmax = (tmax_x < tmax_y) ? tmax_x : tmax_y;

    // Z
    const float tmin_z = (min_z - rayOrg.z) * rayInvDir.z;
    const float tmax_z = (max_z - rayOrg.z) * rayInvDir.z;

    tmin = (tmin > tmin_z) ? tmin : tmin_z;
    tmax = (tmax < tmax_z) ? tmax : tmax_z;

    //
    // Hit include (tmin == tmax) edge case(hit 2D plane).
    //
    if ((tmax > 0.0f) && (tmin <= tmax) && (tmin <= maxT)) {
        *tminOut = tmin;
        *tmaxOut = tmax;
        return true;
    }
    return false; // no hit
}

bool TriangleIsect(float *tInOut, float *uOut, float *vOut,
                   float4 p0, float4 p1, float4 p2,
                   float4 rayOrg, float4 rayDir)
{
    float4 e1, e2;
    float4 p, s, q;

    e1 = p1 - p0;
    e2 = p2 - p0;
    p = cross(rayDir, e2);

    e1.w = 0;
    p.w = 0;

    float invDet;
    float det = dot(e1, p);
    invDet = 1.0f / det;

    s = rayOrg - p0;
    q = cross(s, e1);

    float u = dot(s, p) * invDet;
    float v = dot(q, rayDir) * invDet;
    float t = dot(e2, q) * invDet;

    if (u < 0.0f || u > 1.0f)
        return false;
    if (v < 0.0f || u + v > 1.0f)
        return false;
    if (t < 0.0f || t > *tInOut)
        return false;

    *tInOut = t;
    *uOut = u;
    *vOut = v;
    return true;
}

bool PatchIsect(struct Intersection *isect,
                __global const float *bezierVerts,
                const struct Ray *ray)
{
#if 0
    float4 p0 = (float4)(bezierVerts[0*3], bezierVerts[0*3+1], bezierVerts[0*3+2], 0);
    float4 p1 = (float4)(bezierVerts[3*3], bezierVerts[3*3+1], bezierVerts[3*3+2], 0);
    float4 p2 = (float4)(bezierVerts[12*3], bezierVerts[12*3+1], bezierVerts[12*3+2], 0);
    float4 p3 = (float4)(bezierVerts[15*3], bezierVerts[15*3+1], bezierVerts[15*3+2], 0);

    float4 rayOrg = ray->org;
    float4 rayDir = ray->dir;

    float u, v;
    if (TriangleIsect(&isect->t, &u, &v,
                      p0, p1, p2, rayOrg, rayDir) ||
        TriangleIsect(&isect->t, &u, &v,
                      p1, p2, p3, rayOrg, rayDir)) {
        isect->normal = (float4)(0, 0, -1, 0);
        return true;
    }
    return false;
#else
    struct BezierPatchIntersection bpi;
    struct BezierPatch bp;
    BpConstruct(&bp, bezierVerts);
    BpiConstruct(&bpi, &bp);
    return BpiTest(&bpi, isect, ray, 0, isect->t);
#endif
}

bool TestLeafNode(struct Intersection *isect, // [inout]
                  const struct BVHNode *node,
                  __global const unsigned int *indices,
                  __global const float *bezierVerts,
                  const struct Ray *ray)
{
    unsigned int numPatches = node->data[0];
    unsigned int offset = node->data[1];

    bool hit = false;
    for (unsigned int i = 0; i < numPatches; i++) {
        int faceIdx = indices[i + offset];

        if (PatchIsect(isect, &bezierVerts[faceIdx*16*3], ray)) {
            // Update isect state
            isect->faceID = faceIdx;
            hit = true;
        }
    }
    return hit;
}


__kernel void traverse(__global struct Ray *rays,
                       __global const struct BVHNode *nodes,
                       __global const unsigned int *indices,
                       __global const float *bezierVerts,
                       __global float *image) {

    int gid = get_global_id(0);
    struct Ray ray = rays[gid];

    int nodeStackIndex = 0;
    int nodeStack[5120];
    nodeStack[0] = 0;

    float hitT = FLT_MAX;
    // // Init isect info as no hit
    struct Intersection isect;
    isect.t = hitT;
    isect.u = 0.0f;
    isect.v = 0.0f;
    isect.faceID = -1;

    int dirSign[3];
    dirSign[0] = ray.dir.x < 0.0f ? 1 : 0;
    dirSign[1] = ray.dir.y < 0.0f ? 1 : 0;
    dirSign[2] = ray.dir.z < 0.0f ? 1 : 0;
    ray.dirSign[0] = dirSign[0];
    ray.dirSign[1] = dirSign[1];
    ray.dirSign[2] = dirSign[2];

    ray.invDir.x = 1.0f / ray.dir.x;
    ray.invDir.y = 1.0f / ray.dir.y;
    ray.invDir.z = 1.0f / ray.dir.z;


    float minT, maxT;
    int count = 0;
    while (nodeStackIndex >= 0) {
        ++count;
        int index = nodeStack[nodeStackIndex];
        struct BVHNode node = nodes[index];

        nodeStackIndex--;

        bool hit = IntersectRayAABB(&minT, &maxT, hitT,
                                    node.bmin, node.bmax, ray.org,
                                    ray.invDir, dirSign);

        if (node.flag == 0) { // branch node
            if (hit) {
                int orderNear = dirSign[node.axis];
                int orderFar = 1 - orderNear;

                // Traverse near first.
                nodeStack[++nodeStackIndex] = node.data[orderFar];
                nodeStack[++nodeStackIndex] = node.data[orderNear];
            }
        } else { // leaf node
            if (hit) {
                if (TestLeafNode(&isect, &node, indices, bezierVerts, &ray)) {
                    hitT = isect.t;
                }
            }
        }
    }

    int id = ray.dirSign[3];

    if (hitT < FLT_MAX) {
        float d = fmax(0.0f, dot(ray.dir, isect.normal));
        image[id*4+0] = d;
        image[id*4+1] = d;
        image[id*4+2] = d;
        image[id*4+3] = 1;
    } else {
        image[id*4+0] = 0.1f;
        image[id*4+1] = 0.1f;
        image[id*4+2] = 0.2f;
        image[id*4+3] = 1.0f;
    }
}




