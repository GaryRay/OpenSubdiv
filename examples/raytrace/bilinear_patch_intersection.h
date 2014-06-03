#ifndef __MALLIE_BILINEAR_PATCH_INTERSECTION_H__
#define __MALLIE_BILINEAR_PATCH_INTERSECTION_H__

#include "common.h"
#include "intersection.h"
#include "bounded_intersection.h"

namespace mallie
{
	
	class bilinear_patch_intersection:public bounded_intersection
	{
	public:
		bilinear_patch_intersection(const vector3& p00, const vector3& p10, const vector3& p01, const vector3& p11);
		bilinear_patch_intersection(const vector3 p[4]);
		~bilinear_patch_intersection();
	public:
		bool test    (                    const Ray& r, float dist)const;
		bool test    (Intersection* info, const Ray& r, float dist)const;
		void finalize(Intersection* info, const Ray& r, float dist)const;
		vector3 min()const;
		vector3 max()const;
	private:
		vector3 p_[4];
	};
	
    bool test_bilinear_patch(float* t, float* u, float* v, const vector3f p[4],               float tmin, float tmax);
    bool test_bilinear_patch(double* t, double* u, double* v, const vector3d p[4],          double tmin, double tmax);
    
	bool test_bilinear_patch(float* t, float* u, float* v, const vector3f p[4], const Ray& r, float tmin, float tmax);
	bool test_bilinear_patch(Intersection* info,           const vector3f p[4], const Ray& r, float tmin, float tmax);
	
	
	

	
}

#endif

