#ifndef __BVH_ACCEL_H__
#define __BVH_ACCEL_H__

#include <vector>

#include "common.h"
#include "mesh.h"
#include "intersection.h"

class BVHNode {
public:
  BVHNode() {};
  ~BVHNode() {};

  real bmin[3];
  real bmax[3];

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

///< BVH build option.
struct BVHBuildOptions {
  bool debugPrint;
  real costTaabb;
  int minLeafPrimitives;
  int maxTreeDepth;
  int binSize;

  // Set default value: Taabb = 0.2
  BVHBuildOptions()
      : costTaabb(0.2), minLeafPrimitives(16), maxTreeDepth(256), binSize(64) {}
};

///< BVH build statistics.
struct BVHBuildStatistics {
  int maxTreeDepth;
  int numLeafNodes;
  int numBranchNodes;

  // Set default value: Taabb = 0.2
  BVHBuildStatistics() : maxTreeDepth(0), numLeafNodes(0), numBranchNodes(0) {}
};

class BVHAccel {
public:
  BVHAccel(float uvMargin=0.1f) : _intersectKernel(ORIGINAL),
                                  _uvMargin(uvMargin), _cropUV(true),
                                  _bezierClip(true), _displaceScale(0), _displaceFreq(100) ,
                                  _epsilon(1e-4), _maxLevel(10), _useTriangle(false), _useRayDiffEpsilon(false){}
  ~BVHAccel() {};

  ///< Build BVH for input mesh.
  bool Build(const Mesh *mesh, const BVHBuildOptions &options);

  ///< Get statistics of built BVH tree. Valid after Build()
  BVHBuildStatistics GetStatistics() const { return stats_; }

  ///< Dump built BVH to the file.
  bool Dump(const char *filename);

  /// Load BVH binary
  bool Load(const char *filename);

  ///< Traverse into BVH along ray and find closest hit point if found
  bool Traverse(Intersection &isect, const Mesh *mesh, Ray &ray);

  enum { ORIGINAL, NEW_FLOAT, NEW_SSE, NEW_DOUBLE, OPENCL } IntersectKernel;
  void SetIntersectKernel(int k) {_intersectKernel = k; }
  void SetUVMargin(float margin) { _uvMargin = margin; }
  void SetCropUV(bool flag) {_cropUV = flag;}
  void SetBezierClip(bool flag) {_bezierClip = flag;}
  void SetDisplacement(float scale, float freq) { _displaceScale = scale; _displaceFreq = freq; }
    bool IsGpuKernel() const { return _intersectKernel == OPENCL; }

  void SetEpsilon(double eps){_epsilon=eps;}
  void SetMaxLevel(int level){_maxLevel=level;}
  void SetUseTriangle(bool flag){_useTriangle=flag;}
  void SetUseRayDiffEpsilon(bool flag){_useRayDiffEpsilon=flag;}
  void SetConservativeTest(bool flag){_conservativeTest=flag;}

  const std::vector<BVHNode> &GetNodes() const { return nodes_; }
  const std::vector<unsigned int> &GetIndices() const { return indices_; }

private:
  ///< Builds BVH tree recursively.
  size_t BuildTree(const Mesh *mesh, unsigned int leftIdx,
                   unsigned int rightIdx, int depth);

  BVHBuildOptions options_;
  std::vector<BVHNode> nodes_;
  std::vector<unsigned int> indices_; // max 4G triangles.
  BVHBuildStatistics stats_;

  int _intersectKernel;
  float _uvMargin;
  bool _cropUV;
  bool _bezierClip;
  float _displaceScale;
  float _displaceFreq;

  double _epsilon;
  int    _maxLevel; 
  bool _useTriangle;
  bool _useRayDiffEpsilon;
  bool _conservativeTest;
};

#endif // __BVH_ACCEL_H__
       // vim:set sw=2 ts=2 expandtab:
