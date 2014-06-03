#ifndef __MALLIE_BEZIER_PATCH_INTERSECTION_H__
#define __MALLIE_BEZIER_PATCH_INTERSECTION_H__

//#include "bounded_intersection.h"

#include "common.h"
#include "intersection.h"
#include "bezier_patch.hpp"
#include "bounded_intersection.h"
#include "vector3.h"

namespace mallie{
	
	class bezier_patch_intersection:public bounded_intersection
	{
	public:
		bezier_patch_intersection(const bezier_patch<vector3>& patch);
		bezier_patch_intersection(const bezier_patch<vector3>& patch, float u0, float u1, float v0, float v1);
		~bezier_patch_intersection();
	public:
		bool test    (                    const Ray& r, float tmin, float tmax)const;
		bool test    (Intersection* info, const Ray& r, float tmin, float tmax)const;
	public:
		bool test    (                    const Ray& r, float dist)const;
		bool test    (Intersection* info, const Ray& r, float dist)const;
		void finalize(Intersection* info, const Ray& r, float dist)const;
		vector3 min()const;
		vector3 max()const;
	public:
		void finalize(float u, float v, float t, Intersection* info, const Ray& r, float dist)const;
    protected:
        bool test_internal    (                 const Ray& r, float tmin, float tmax)const;
        bool test_internal    (Intersection* info, const Ray& r, float tmin, float tmax)const;
	protected:
		bezier_patch<vector3> patch_;
		float urange_[2];
		float vrange_[2];
		vector3 min_;
		vector3 max_;
		float eps_;
	};
	
}

#endif

