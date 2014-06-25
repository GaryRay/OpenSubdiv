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
    float org[4];
    float dir[4];
    float invDir[4];
    int dirSign[4];
};

struct BezierPatch {
    float4 cp[16];
};

struct BezierPatchIntersection {
    struct BezierPatch patch;
    float4 min;
    float4 max;
    float eps;
};

struct RangeAABB {
    float tmin, tmax;
};


// --------------------------------------------------------------------
// prototypes
bool IntersectRayAABB(float *tminOut, float *tmaxOut,
                      float maxT, float *bmin, float *bmax,
                      float *rayOrg, float *rayInvDir, int *rayDirSign);
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
void BpiConstruct(struct BezierPatchIntersection *isect,
                  const struct BezierPatch *patch);
bool BpiTest(struct BezierPatchIntersection *isect,
             struct Intersection *info,
             const struct Ray *r,
             float tmin, float tmax);
bool IntersectAABB(struct RangeAABB *rng, float4 min, float4 max,
                   const struct Ray *r, float tmin, float tmax);
// --------------------------------------------------------------------
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
    *max = patch->cp[1];
#pragma unroll
    for (int i = 1; i < 16; ++i) {
        *min = fmin(*min, patch->cp[i]);
        *max = fmax(*max, patch->cp[i]);
    }
    *min -= (float4)(eps, eps, eps, 0);
    *max += (float4)(eps, eps, eps, 0);
}

void BpiConstruct(struct BezierPatchIntersection *isect,
                  const struct BezierPatch *patch)
{
    BpGetMinMax(patch, &isect->min, &isect->max, 0.01f);
    isect->eps = 1.0e-4f;
}

bool IntersectAABB(struct RangeAABB *rng, float4 min, float4 max,
                   const struct Ray *r, float tmin, float tmax)
{
    float4 box[2];
    box[0] = min;
    box[1] = max;
    int4 sign = (int4)(r->dirSign[0], r->dirSign[1], r->dirSign[2], 0);

    float4 org = (float4)(r->org[0], r->org[1], r->org[2], 0);
    float4 idir = (float4)(r->invDir[0], r->invDir[1], r->invDir[2], 0);

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
        info->t = 0;

        // TODO.
        return true;
    }
    return false;
}


// --------------------------------------------------------------------

bool IntersectRayAABB(float *tminOut, float *tmaxOut,
                      float maxT, float *bmin, float *bmax,
                      float *rayOrg, float *rayInvDir, int *rayDirSign)
{
    float tmin, tmax;

    const float min_x = rayDirSign[0] * bmax[0] + (1-rayDirSign[0]) * bmin[0];
    const float min_y = rayDirSign[1] * bmax[1] + (1-rayDirSign[1]) * bmin[1];
    const float min_z = rayDirSign[2] * bmax[2] + (1-rayDirSign[2]) * bmin[2];
    const float max_x = rayDirSign[0] * bmin[0] + (1-rayDirSign[0]) * bmax[0];
    const float max_y = rayDirSign[1] * bmin[1] + (1-rayDirSign[1]) * bmax[1];
    const float max_z = rayDirSign[2] * bmin[2] + (1-rayDirSign[2]) * bmax[2];

    // X
    const float tmin_x = (min_x - rayOrg[0]) * rayInvDir[0];
    const float tmax_x = (max_x - rayOrg[0]) * rayInvDir[0];

    // Y
    const float tmin_y = (min_y - rayOrg[1]) * rayInvDir[1];
    const float tmax_y = (max_y - rayOrg[1]) * rayInvDir[1];

    tmin = (tmin_x > tmin_y) ? tmin_x : tmin_y;
    tmax = (tmax_x < tmax_y) ? tmax_x : tmax_y;

    // Z
    const float tmin_z = (min_z - rayOrg[2]) * rayInvDir[2];
    const float tmax_z = (max_z - rayOrg[2]) * rayInvDir[2];

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

    float4 rayOrg = (float4)(ray->org[0], ray->org[1], ray->org[2], 0);
    float4 rayDir = (float4)(ray->dir[0], ray->dir[1], ray->dir[2], 0);

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
                       int stepIndex,
                       int step,
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
    dirSign[0] = ray.dir[0] < 0.0f ? 1 : 0;
    dirSign[1] = ray.dir[1] < 0.0f ? 1 : 0;
    dirSign[2] = ray.dir[2] < 0.0f ? 1 : 0;

    ray.dirSign[0] = dirSign[0];
    ray.dirSign[1] = dirSign[1];
    ray.dirSign[2] = dirSign[2];

    float rayInvDir[3];
    rayInvDir[0] = 1.0f / ray.dir[0];
    rayInvDir[1] = 1.0f / ray.dir[1];
    rayInvDir[2] = 1.0f / ray.dir[2];

    ray.invDir[0] = rayInvDir[0];
    ray.invDir[1] = rayInvDir[1];
    ray.invDir[2] = rayInvDir[2];


    float minT, maxT;
    int count = 0;
    while (nodeStackIndex >= 0) {
        ++count;
        int index = nodeStack[nodeStackIndex];
        struct BVHNode node = nodes[index];

        nodeStackIndex--;

        bool hit = IntersectRayAABB(&minT, &maxT, hitT,
                                    node.bmin, node.bmax, ray.org,
                                    rayInvDir, dirSign);

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
        image[id*4+0] = count;
        image[id*4+1] = count * 0.5f;
        image[id*4+2] = count * 0.2f;
        image[id*4+3] = 1;
    } else {
        image[id*4+0] = 0;
        image[id*4+1] = 0;
        image[id*4+2] = 0;
        image[id*4+3] = 1;
    }
}




