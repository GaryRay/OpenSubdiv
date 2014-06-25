#include <cstdio>
#include <osd/opencl.h>
#include "clTracer.h"

static const char *clSource =
#include "clKernel.gen.h"
    ;

#define CL_CHECK_ERROR(x, ...) {      \
        if (x != CL_SUCCESS)          \
        { printf("ERROR %d : ", x);   \
            printf(__VA_ARGS__);} }

static cl_kernel buildKernel(cl_program prog, const char * name) {

    cl_int ciErr;
    cl_kernel k = clCreateKernel(prog, name, &ciErr);

    if (ciErr != CL_SUCCESS) {
        printf("Can't compile kernel\n");
    }
    return k;
}

static bool initCL(cl_context *clContext, cl_command_queue *clQueue)
{
    cl_platform_id cpPlatform = 0;
    cl_uint num_platforms;
    cl_int ciErrNum = clGetPlatformIDs(0, NULL, &num_platforms);
    if (ciErrNum != CL_SUCCESS) {
        printf("Error %d in clGetPlatformIDs call.\n", ciErrNum);
        return false;
    }
    if (num_platforms == 0) {
        printf("No OpenCL platform found.\n");
        return false;
    }
    cl_platform_id *clPlatformIDs = new cl_platform_id[num_platforms];
    ciErrNum = clGetPlatformIDs(num_platforms, clPlatformIDs, NULL);
    char chBuffer[1024];
    for (cl_uint i = 0; i < num_platforms; ++i) {
        ciErrNum = clGetPlatformInfo(clPlatformIDs[i], CL_PLATFORM_NAME,
                                     1024, chBuffer,NULL);
        if (ciErrNum == CL_SUCCESS) {
            cpPlatform = clPlatformIDs[i];
        }
    }

    cl_uint numDevices = 0;
    clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, 0, NULL, &numDevices);
    if (numDevices == 0) {
        printf("No sharable devices.\n");
        return false;
    }
    cl_device_id *clDevices = new cl_device_id[numDevices];
    clGetDeviceIDs(cpPlatform, CL_DEVICE_TYPE_GPU, numDevices, clDevices, NULL);

    // cl_context_properties props[] = {
    //     CL_CONTEXT_PLATFORM, (cl_context_properties)cpPlatform,
    //     0
    // };

    *clContext = clCreateContext(NULL, numDevices, clDevices, NULL, NULL, &ciErrNum);
    if (ciErrNum != CL_SUCCESS) {
        printf("Error %d in clCreateContext\n", ciErrNum);
        delete[] clDevices;
        return false;
    }

    *clQueue = clCreateCommandQueue(*clContext, clDevices[0], 0, &ciErrNum);
    delete[] clDevices;
    if (ciErrNum != CL_SUCCESS) {
        printf("Error %d in clCreateCommandQueue\n", ciErrNum);
        return false;
    }
    return true;
}

CLTracer::CLTracer() :
    _rays(0), _bvhNodes(0), _bvhIndices(0), _bezierVerts(0)
{
    compile();
}

void
CLTracer::compile()
{
    initCL(&_clContext, &_clQueue);

    cl_int ciErrNum;
    const char *sources[] = { clSource };

    _clProgram = clCreateProgramWithSource(_clContext, 1, sources, 0, &ciErrNum);
    CL_CHECK_ERROR(ciErrNum, "clCreateProgramWithSource\n");

    ciErrNum = clBuildProgram(_clProgram, 0, NULL, NULL, NULL, NULL);
    CL_CHECK_ERROR(ciErrNum, "clBuildProgram\n");

    if (ciErrNum != CL_SUCCESS) {
        cl_int numDevices = 0;
        clGetContextInfo(_clContext, CL_CONTEXT_NUM_DEVICES,
                         sizeof(cl_uint), &numDevices, NULL);
        cl_device_id *devices = new cl_device_id[numDevices];
        clGetContextInfo(_clContext, CL_CONTEXT_DEVICES,
                         sizeof(cl_device_id)*numDevices, devices, NULL);
        for (int i = 0; i < numDevices; ++i) {
            char cBuildLog[10240];
            clGetProgramBuildInfo(_clProgram, devices[i], CL_PROGRAM_BUILD_LOG,
                                  sizeof(cBuildLog), cBuildLog, NULL);
            printf(cBuildLog);
        }
        delete[] devices;
        return;
    }

    _kernel = buildKernel(_clProgram, "traverse");
}

void
CLTracer::SetBVH(BVHAccel const &accel)
{
    std::vector<BVHNode> const &node = accel.GetNodes();
    std::vector<unsigned int> const &indices= accel.GetIndices();

    if (_bvhNodes) clReleaseMemObject(_bvhNodes);
    if (_bvhIndices) clReleaseMemObject(_bvhIndices);

    cl_int ciErrNum;
    _bvhNodes = clCreateBuffer(
        _clContext,
        CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
        node.size()*sizeof(BVHNode),
        const_cast<BVHNode*>(&node[0]), &ciErrNum);
    CL_CHECK_ERROR(ciErrNum, "clCreateBuffer\n");

    _bvhIndices = clCreateBuffer(
        _clContext,
        CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
        indices.size()*sizeof(unsigned int),
        const_cast<unsigned int*>(&indices[0]), &ciErrNum);
    CL_CHECK_ERROR(ciErrNum, "clCreateBuffer\n");
}

void
CLTracer::SetBezierVertices(const float *bezierVerts, int size)
{
    if (_bezierVerts) clReleaseMemObject(_bezierVerts);

    cl_int ciErrNum;
    _bezierVerts = clCreateBuffer(
        _clContext,
        CL_MEM_READ_ONLY|CL_MEM_COPY_HOST_PTR,
        size,
        const_cast<float*>(bezierVerts), &ciErrNum);
    CL_CHECK_ERROR(ciErrNum, "clCreateBuffer\n");
}

void
CLTracer::Traverse(int width, int height, const CLRay *rays, float *image)
{
    cl_int ciErrNum;

    if (_rays) clReleaseMemObject(_rays);
    _rays = clCreateBuffer(
        _clContext,
        CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
        width*height*sizeof(CLRay),
        const_cast<CLRay*>(rays), &ciErrNum);
    CL_CHECK_ERROR(ciErrNum, "clCreateBuffer\n");

    cl_mem dst = clCreateBuffer(
        _clContext,
        CL_MEM_READ_WRITE|CL_MEM_COPY_HOST_PTR,
        width*height*4*sizeof(float),
        image, &ciErrNum);
    CL_CHECK_ERROR(ciErrNum, "clCreateBuffer\n");

    clSetKernelArg(_kernel, 0, sizeof(cl_mem), &_rays);
    clSetKernelArg(_kernel, 1, sizeof(cl_mem), &_bvhNodes);
    clSetKernelArg(_kernel, 2, sizeof(cl_mem), &_bvhIndices);
    clSetKernelArg(_kernel, 3, sizeof(cl_mem), &_bezierVerts);
    clSetKernelArg(_kernel, 4, sizeof(cl_mem), &dst);

    size_t globalWorkSize[1] = { (size_t)width*height };
    ciErrNum = clEnqueueNDRangeKernel(_clQueue,
                                      _kernel, 1, NULL, globalWorkSize,
                                      NULL, 0, NULL, NULL);
    CL_CHECK_ERROR(ciErrNum, "clEnqueueNDRangeKernel\n");

    ciErrNum = clEnqueueReadBuffer(_clQueue, dst, true,
                                   0, width*height*sizeof(float)*4, image,
                                   0, NULL, NULL);
    CL_CHECK_ERROR(ciErrNum, "clEnqueueReadBuffer\n");

    clReleaseMemObject(dst);
}
