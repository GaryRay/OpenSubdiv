#include "bilinear_patch_intersection.h"
#include "intersection.h"

#include "vector3.h"

#include <math.h>
#include <string.h>
#include <assert.h>

#include <utility>
#include <algorithm>
#include <limits>

namespace mallie{

    static const float EPS = std::numeric_limits<float>::epsilon();
    
    struct range_AABB{
        float tmin;
        float tmax;
    };

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
    
    static inline bool intersectAABB(range_AABB* rng, const vector3 & min, const vector3 & max, const Ray & r, float tmin, float tmax)
    {
        //int phase = r.dirSign;
        int sign[3];// = {(phase>>0)&1,(phase>>1)&1,(phase>>2)&1};
        memcpy(sign, r.dirSign, sizeof(int)*3);
        vector3 box[2] = {min,max};
        vector3 org  = Conv(r.org);
        vector3 idir = Conv(r.invDir);

        for(int i = 0;i<3;i++){
            tmin = std::max<float>(tmin,(box[  sign[i]][i]-org[i])*idir[i]);
            tmax = std::min<float>(tmax,(box[1-sign[i]][i]-org[i])*idir[i]);
        }
        rng->tmin = tmin;
        rng->tmax = tmax;
		return rng->tmin<=rng->tmax;

    }
    
	
	bilinear_patch_intersection::bilinear_patch_intersection(const vector3& p00, const vector3& p10, const vector3& p01, const vector3& p11)
	{
		p_[0] = p00;
		p_[1] = p10;
		p_[2] = p01;
		p_[3] = p11;
	}
	
	bilinear_patch_intersection::bilinear_patch_intersection(const vector3 p[4])
	{
		memcpy(p_, p, sizeof(vector3)*4);
	}
	
	bilinear_patch_intersection::~bilinear_patch_intersection()
	{
		;//
	}
	
	bool bilinear_patch_intersection::test    (                 const Ray& r, float dist)const
	{
		return test_bilinear_patch(NULL, p_, r, 0, dist);
	}
	
	bool bilinear_patch_intersection::test    (Intersection* info, const Ray& r, float dist)const
	{
		if(test_bilinear_patch(info, p_, r, 0, dist)){
			//info->p_intersection = this;
			finalize(info, r, info->t);
			return true;
		}
		return false;
	}

    static
    inline vector3 leap(const vector3& a, const vector3& b, float t)
    {
        return a*(1.0-t) + b*t;
    }

	
	void bilinear_patch_intersection::finalize(Intersection* info, const Ray& r, float dist)const
	{
		real u = info->u;
		real v = info->v;
		
		vector3 U = normalize(leap(p_[1],p_[3],v)-leap(p_[0],p_[2],v));
		vector3 V = normalize(leap(p_[2],p_[3],u)-leap(p_[0],p_[1],u));
		
		vector3 N = cross(U, V);
		
		info->position = r.org + dist * r.dir;
		info->normal = Conv(N);
		info->geometricNormal = Conv(N);
		info->tangent = Conv(U);
		info->binormal = Conv(V);
	}
	
	vector3 bilinear_patch_intersection::min()const
	{
		vector3 m = p_[0]; 
		for(int i=1;i<4;i++)
		{
			for(int j=0;j<3;j++){
				if(m[j]>p_[i][j])m[j]=p_[i][j];
			}
		}
        return m - vector3(EPS, EPS, EPS);
	}
	
	vector3 bilinear_patch_intersection::max()const
	{
		vector3 m = p_[0]; 
		for(int i=1;i<4;i++)
		{
			for(int j=0;j<3;j++){
				if(m[j]<p_[i][j])m[j]=p_[i][j];
			}
		}
		return m + vector3(EPS, EPS, EPS);
	}
	
    template<class FLOAT>
	static
	bool solve_bilinear_patch(FLOAT* t, FLOAT* u, FLOAT* v, FLOAT tmin, FLOAT tmax, FLOAT tt, FLOAT uu, FLOAT vv)
	{
        if(tmin<=tt&&tt<=tmax){
            *t = tt;
            *u = uu;
            *v = vv;
            return true;
        }
		return false;
	}

	template<class FLOAT>
    static
    int solve2(FLOAT root[2], const FLOAT coeff[3])
    {
        FLOAT A = coeff[0];
        FLOAT B = coeff[1];
        FLOAT C = coeff[2];
        
        if(fabs(A)<=std::numeric_limits<FLOAT>::min()){
        	if(fabs(B)<=std::numeric_limits<FLOAT>::min())return 0;
            FLOAT x = -C/B;
            root[0] = x;
            return 1;
        }else{  
            FLOAT D = B*B-4*A*C;
            if(D<0){
                return 0;
            }else if(D==0){
                FLOAT x = -0.5*B/A;
                root[0] = x;
                return 1;
            }else{
                FLOAT x1 = (fabs(B) + sqrt(D)) / (2.0 * A);
                if ( B >= 0.0 ) {
                    x1 = -x1;
                }
                FLOAT x2 = C / (A * x1);
                if(x1>x2)std::swap<FLOAT>(x1,x2);
                root[0] = x1;
                root[1] = x2;
                return 2;
            }
        }
    }
    
    template<class FLOAT>
    static
    inline FLOAT compute_u(FLOAT A1, FLOAT A2, FLOAT B1, FLOAT B2, FLOAT C1, FLOAT C2, FLOAT D1, FLOAT D2, FLOAT v)
    {
        //return div((v*(C1-C2)+(D1-D2)),(v*(A2-A1)+(B2-B1)));
        FLOAT a = v*A2+B2;
        FLOAT b = v*(A2-A1)+B2-B1;
        if(fabs(b)>=fabs(a)){
            return (v*(C1-C2)+D1-D2)/b;
        }else{
            return (-v*C2-D2)/a;
        }
    }
    
    template<class FLOAT>
    static
    inline FLOAT compute_t(FLOAT a, FLOAT b, FLOAT c, FLOAT d, FLOAT iq, FLOAT u, FLOAT v)
    {
        return ((u*v)*a + u*b + v*c + d)*iq;
    }
    
    template<class V>
	static
	void get_minmax(V& min, V& max, const V p[4])
	{
        min = max = p[0];
        for(int i = 1;i<4;i++)
        {
            for(int j=0;j<3;j++)
            {
                if(min[j]>p[i][j])min[j]=p[i][j];
                if(max[j]<p[i][j])max[j]=p[i][j];
            }
        }
	}

	
#define CALC_U(v) compute_u(A1,A2,B1,B2,C1,C2,D1,D2, v)
#define CALC_T(u,v) compute_t(a[nPlane],b[nPlane],c[nPlane],d[nPlane], iq[nPlane], u, v)
	
	bool test_bilinear_patch(float* t, float* u, float* v, const vector3f p[4], const Ray& r, float tmin, float tmax)
	{
		const vector3& p00 = p[0];
		const vector3& p10 = p[1];
		const vector3& p01 = p[2];
		const vector3& p11 = p[3];
		
		vector3 q  = Conv(r.dir);
        vector3 iq = Conv(r.invDir);
        
        int nPlane = 0;
        {
            float DD[3];
            DD[0] = fabs(iq[0]);
            DD[1] = fabs(iq[1]);
            DD[2] = fabs(iq[2]);
        
            if(DD[nPlane]>DD[1])nPlane=1;
            if(DD[nPlane]>DD[2])nPlane=2;
        }
		
		vector3 a = p11-p10-p01+p00;
		vector3 b = p10-p00;
		vector3 c = p01-p00;
		vector3 d = p00 - Conv(r.org);
		
        //xz-zx
		float A1 = a[0]*q[2]-a[2]*q[0];
		float B1 = b[0]*q[2]-b[2]*q[0];
		float C1 = c[0]*q[2]-c[2]*q[0];
		float D1 = d[0]*q[2]-d[2]*q[0];
		
        //yz-zy
		float A2 = a[1]*q[2]-a[2]*q[1];
		float B2 = b[1]*q[2]-b[2]*q[1];
		float C2 = c[1]*q[2]-c[2]*q[1];
		float D2 = d[1]*q[2]-d[2]*q[1];
		
		float F1 = A2*C1-A1*C2;
		float F2 = A2*D1-A1*D2+B2*C1-B1*C2;
		float F3 = B2*D1-B1*D2;

        float coeff[] = {F1,F2,F3};

        float root[2] = {};
        int nRet = solve2(root,coeff); 
	
        if(nRet)
        {
            bool bRet = false;
            for(int i=0;i<nRet;i++)
            {
                float vv = root[i];
			    if(0<=vv && vv<=1){
				    float uu = CALC_U(vv);
				    if(0<=uu && uu<=1){
                        float tt = CALC_T(uu,vv);
					    if(solve_bilinear_patch(t, u, v, tmin, tmax, tt, uu, vv)){
						    tmax = *t; 
						    bRet = true;
					    }
				    }
			    }
            }
            return bRet;
        }
        return false;
	}

#undef CALC_T
#undef CALC_U
    
#define CALC_U(v) compute_u(A1,A2,B1,B2,C1,C2,D1,D2, v)
#define CALC_T(u,v) compute_t(a[nPlane],b[nPlane],c[nPlane],d[nPlane], (FLOAT)1.0, u, v)
	
	template<class FLOAT, class V>
    bool test_bilinear_patch_(FLOAT* t, FLOAT* u, FLOAT* v, const V p[4], FLOAT tmin, FLOAT tmax)
	{
		const V& p00 = p[0];
		const V& p10 = p[1];
		const V& p01 = p[2];
		const V& p11 = p[3];

		/*
		V min;
		V max;
		get_minmax(min, max, p);
		if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<tmin || tmax<min[2])return false;//z
		*/


        static const int nPlane = 2;
		
		V a = p11-p10-p01+p00;
		V b = p10-p00;
		V c = p01-p00;
		V d = p00;
		
        //xz-zx
		FLOAT A1 = a[0];
		FLOAT B1 = b[0];
		FLOAT C1 = c[0];
		FLOAT D1 = d[0];
		
        //yz-zy
		FLOAT A2 = a[1];
		FLOAT B2 = b[1];
		FLOAT C2 = c[1];
		FLOAT D2 = d[1];
		
		FLOAT F1 = A2*C1-A1*C2;
		FLOAT F2 = A2*D1-A1*D2+B2*C1-B1*C2;
		FLOAT F3 = B2*D1-B1*D2;
        
        FLOAT coeff[] = {F1,F2,F3};
        
        FLOAT root[2] = {};
        int nRet = solve2(root,coeff);
        
        if(nRet)
        {
            bool bRet = false;
            for(int i=0;i<nRet;i++)
            {
                FLOAT vv = root[i];
			    if(0<=vv && vv<=1){
				    FLOAT uu = CALC_U(vv);
				    if(0<=uu && uu<=1){
                        FLOAT tt = CALC_T(uu,vv);
					    if(solve_bilinear_patch(t, u, v, tmin, tmax, tt, uu, vv)){
						    tmax = *t; 
						    bRet = true;
					    }
				    }
			    }
            }
            return bRet;
        }
        return false;
        
    }
    
#undef CALC_T
#undef CALC_U

	bool test_bilinear_patch(float* t, float* u, float* v, const vector3f p[4], float tmin, float tmax){
		return test_bilinear_patch_(t, u, v, p, tmin, tmax);
	}
	bool test_bilinear_patch(double* t, double* u, double* v, const vector3d p[4], double tmin, double tmax){
		return test_bilinear_patch_(t, u, v, p, tmin, tmax);
	}

	
	bool test_bilinear_patch(Intersection* info, const vector3 p[4], const Ray& r, float tmin, float tmax)
	{
		float t = 0;
		float u = 0;
		float v = 0;
		if(test_bilinear_patch(&t,&u,&v, p, r, tmin, tmax))
		{
			if(info)
			{
				info->t = t;
				info->texcoord[0] = u;
				info->texcoord[1] = v;
			}
			return true;
		}
		return false;
	}
	
}

