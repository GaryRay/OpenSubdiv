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
    int level;
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
void BpRotateU(struct BezierPatch *patch);

void BpiConstruct(struct BezierPatchIntersection *isect,
                  const struct BezierPatch *patch);
bool BpiTest(struct BezierPatchIntersection *isect,
             struct Intersection *info,
             const struct Ray *r,
             float tmin, float tmax);
bool IntersectAABB(struct RangeAABB *rng, float4 min, float4 max,
                   const struct Ray *r, float tmin, float tmax);
bool BpiTestInternal(struct BezierPatchIntersection *isect,
                     struct Intersection *info, const struct Ray *r,
                     float tmin, float tmax);
bool BpiTestBezierPatch(struct BezierPatchIntersection *this,
                        struct UVT *uvt, struct BezierPatch *patch,
                        float tmin, float tmax, float eps);
bool BpiTestBezierClipU(struct BezierPatchIntersection *this,
                        struct UVT *info, struct BezierPatch *patch,
                        float u0, float u1, float v0, float v1,
                        float zmin, float zmax, int level, int max_level, float eps);

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
    struct BezierPatch patch = isect->patch;
    BpTransform(&patch, &mat);
    if (BpiTestBezierPatch(isect, &uvt, &patch, tmin, tmax, isect->eps)) {
        float t = uvt.t;
        float u = uvt.u;
        float v = uvt.v;

        //        u = isect->uRange[0]*(1-u) + isect->uRange[1]*u;//global
        //        v = isect->vRange[0]*(1-v) + isect->vRange[1]*v;//global
        info->t = t;
        info->u = u;
        info->v = v;
        //        info->clipLevel = uvt.level;
#if 0
        {
            float4 du = _patch.EvaluateDu(u,v);
            ValueType dv = _patch.EvaluateDv(u,v);
            ValueType normal = cross(du,dv);
            normal.normalize();
            info->normal = real3(normal[0], normal[1], normal[2]);
            info->geometricNormal = real3(normal[0], normal[1], normal[2]);
            //                info->tangent  = Conv(U);
            //                info->binormal = Conv(V);
        }
#endif
        return true;
    }
    return false;
}

void BpRotateU(struct BezierPatch *patch)
{
    float4 dx = patch->cp[12] - patch->cp[0] + patch->cp[15] - patch->cp[3];
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

bool BpiTestBezierClipU(struct BezierPatchIntersection *this,
                        struct UVT *info, struct BezierPatch *patch,
                        float u0, float u1, float v0, float v1,
                        float zmin, float zmax, int level, int max_level, float eps)
{
    struct BezierPatch tpatch = *patch;
    BpRotateU(&tpatch);

    float4 min, max;
    BpGetMinMax(&tpatch, &min, &max, eps*1e-3f);

    if (0 < min[0] || max[0] < 0) return false;//x
    if (0 < min[1] || max[1] < 0) return false;//y
    if (max[2] < zmin || zmax < min[2]) return false;//z

    // TODO
    info->t = 0;
    return true;
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
    int maxLevel = 10;

    return BpiTestBezierClipU(this, info, patch, 0, 1, 0, 1, zmin, zmax, 0, maxLevel, eps);
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
        image[id*4+0] = isect.u;
        image[id*4+1] = isect.v;
        image[id*4+2] = 1;
        image[id*4+3] = 1;
    } else {
        image[id*4+0] = 0;
        image[id*4+1] = 0;
        image[id*4+2] = 0;
        image[id*4+3] = 1;
    }
}




