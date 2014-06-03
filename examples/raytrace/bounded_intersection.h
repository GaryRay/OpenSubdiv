#ifndef __MALLIE_BOUNDED_INTERSECTION_H__
#define __MALLIE_BOUNDED_INTERSECTION_H__

#include "vector3.h"
#include "intersection.h"
#include "common.h"

namespace mallie{

	class bounded_intersection{
	public:
		virtual ~bounded_intersection(){}
		virtual bool test(const Ray & r, float dist)const=0;
		virtual bool test(Intersection * info, const Ray & r, float dist)const=0;
		virtual void finalize(Intersection * info, const Ray & r, float dist)const=0;
		
		virtual vector3 min()const=0;
		virtual vector3 max()const=0;
	};
	
	typedef bounded_intersection intersector;//

}

#endif
