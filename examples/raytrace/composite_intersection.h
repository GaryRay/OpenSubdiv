#ifndef __COMPOSITE_INTERSECTION_H__
#define __COMPOSITE_INTERSECTION_H__

#include "intersection.h"
#include "bounded_intersection.h"
#include "count_ptr.hpp"

namespace mallie{

	class composite_intersection:public bounded_intersection{
	public:
		virtual bool test(const Ray& r, float dist)const = 0;
		virtual bool test(Intersection* info, const Ray& r, float dist)const = 0;
		virtual void finalize(Intersection* info, const Ray& r, float dist)const = 0;// -- //NO IMPLEMENT!
		
		virtual vector3 min()const = 0;
		virtual vector3 max()const = 0;
	public:
		virtual void construct() = 0;
		virtual void add(const auto_count_ptr<intersector>& inter) = 0;
	public:
		virtual ~composite_intersection(){}
	};
	
}

#endif
