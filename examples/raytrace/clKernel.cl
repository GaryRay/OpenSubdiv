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
}

bool TestLeafNode(struct Intersection *isect, // [inout]
                  const struct BVHNode *node,
                  __global const unsigned int *indices,
                  __global const float *bezierVerts,
                  const struct Ray *ray)
{
    unsigned int numTriangles = node->data[0];
    unsigned int offset = node->data[1];

    float t = isect->t;

    bool hit = false;
    for (unsigned int i = 0; i < numTriangles; i++) {
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
    isect.u = 0.0;
    isect.v = 0.0;
    isect.faceID = -1;

    int dirSign[3];
    dirSign[0] = ray.dir[0] < 0.0f ? 1 : 0;
    dirSign[1] = ray.dir[1] < 0.0f ? 1 : 0;
    dirSign[2] = ray.dir[2] < 0.0f ? 1 : 0;

    float rayInvDir[3];
    rayInvDir[0] = 1.0f / ray.dir[0];
    rayInvDir[1] = 1.0f / ray.dir[1];
    rayInvDir[2] = 1.0f / ray.dir[2];

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
    if (hitT < FLT_MAX) {
        image[gid*4+0] = count;
        image[gid*4+1] = count * 0.5;
        image[gid*4+2] = count * 0.2;
        image[gid*4+3] = 1;
    } else {
        image[gid*4+0] = 0;
        image[gid*4+1] = 0;
        image[gid*4+2] = 0;
        image[gid*4+3] = 1;
    }
}




