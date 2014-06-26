#ifndef CL_TRACER_H
#define CL_TRACER_H

#include <osd/opencl.h>
#include "bvh_accel.h"

struct CLRay {
    CLRay() { };
    CLRay(Ray const &ray, int id) {
        for (int i = 0; i < 3; ++i) {
            org[i] = ray.org[i];
            dir[i] = ray.dir[i];
            invDir[i] = ray.invDir[i];
            dirSign[i] = ray.dirSign[i];
        }
        dirSign[3] = id;
    }
    float org[4];
    float dir[4];
    float invDir[4];
    int dirSign[4];
};

class CLTracer {
public:
    CLTracer();
    ~CLTracer();

    void SetBVH(BVHAccel const &accel);
    void SetBezierVertices(const float *bezierVerts, int size);
    void SetImageSize(int width, int height);

    void Traverse(const CLRay *rays, int step, float *image);

private:
    void compile();

    cl_context _clContext;
    cl_command_queue _clQueue;
    cl_program _clProgram;
    cl_kernel _kernel;

    cl_mem _rays;
    cl_mem _image;
    cl_mem _bvhNodes;
    cl_mem _bvhIndices;
    cl_mem _bezierVerts;

    int _width;
    int _height;
};

#endif  // CL_TRACER_H
