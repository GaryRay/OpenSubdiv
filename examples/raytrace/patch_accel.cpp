#include "patch_accel.h"

#include "bezier/bezier_patch.hpp"
#include "bezier/bezier_patch_intersection.h"
#include "bezier/bvh_composite_intersection.h"


namespace mallie
{
	static
    inline vector3 Conv(const real3& r)
    {
        return vector3(r.x, r.y, r.z); 
    }

    static
    inline real3 Conv(const vector3& r)
    {
        return real3(r[0], r[1], r[2]); 
    }


	class PatchAccelImp {
	public:
		PatchAccelImp(std::vector< bezier_patch<vector3> >& v)
		{
			bvh_ = new bvh_composite_intersection();
			size_t sz = v.size();
			for(size_t i=0;i<sz;i++)
			{
				bvh_->add(new bezier_patch_intersection(v[i]));
			}
			bvh_->construct();
		}
		~PatchAccelImp()
		{
			if(bvh_)delete bvh_;
		}
		bool Traverse(Intersection &isect, const Ray &ray)
		{
			return bvh_->test(&isect, ray, 1e+8);
		}
		void GetBoundingBox(real3 &bmin, real3 &bmax)
		{
			vector3 min = bvh_->min();
			vector3 max = bvh_->max();

			bmin = Conv(min);
			bmax = Conv(max);
		}
	protected:
		bvh_composite_intersection* bvh_;
	};

	PatchAccel::PatchAccel(std::vector< bezier_patch<vector3> >& v)
	{
		imp_ = new PatchAccelImp(v);
	}
	PatchAccel::~PatchAccel()
	{
		delete imp_;
	}
	bool PatchAccel::Traverse(Intersection &isect, const Ray &ray)
	{
		imp_->Traverse(isect, ray);
	}

	void PatchAccel::GetBoundingBox(real3 &bmin, real3 &bmax)
	{
		imp_->GetBoundingBox(bmin,bmax);
	}


}
