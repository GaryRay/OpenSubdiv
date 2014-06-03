#ifndef __BVH_COMPOSITE_INTERSECTION_H__
#define __BVH_COMPOSITE_INTERSECTION_H__

#include "composite_intersection.h"
#include "count_ptr.hpp"

namespace mallie{

	class bvh_composite_intersection_imp;

	class bvh_composite_intersection:public composite_intersection{
	public:
		bvh_composite_intersection();
		~bvh_composite_intersection();
	public:
		bool test(const Ray& r, float dist)const;
		bool test(Intersection* info, const Ray& r, float dist)const;
		void finalize(Intersection* info, const Ray& r, float dist)const;// -- //NO IMPLEMENT!
		
		vector3 min()const;
		vector3 max()const;
	public:
		void construct();
		void add(const auto_count_ptr<intersector>& inter);
	private:
		bvh_composite_intersection_imp* imp_;        
	};

}

#endif