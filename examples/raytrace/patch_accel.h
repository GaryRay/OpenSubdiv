#ifndef __PATCH_ACCEL_H__
#define __PATCH_ACCEL_H__

#include "common.h"
#include "vector3.h"
#include "bezier_patch.hpp"
#include "intersection.h"

namespace mallie
{

class PatchAccelImp;
class PatchAccel {
public:
	PatchAccel(std::vector< bezier_patch<vector3> >& v);
	~PatchAccel();
	bool Traverse(Intersection &isect, const Ray &ray);
	void GetBoundingBox(real3 &bmin, real3 &bmax);
protected:
	PatchAccelImp* imp_;
};

}

#endif