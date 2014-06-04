#include "bezier_patch_intersection.h"
#include "bilinear_patch_intersection.h"

#include "vector2.h"
#include "vector3.h"
#include "matrix3.h"
#include "matrix4.h"
#include "common.h"

#include <iostream>
#include <utility>
#include <limits>

//
#define DIRECT_BILINEAR 1
#define USE_BEZIERCLIP 1
#define USE_COARSESORT 1

#define USE_FLOAT 1
#define USE_SSE 0

#if USE_SSE
#include "sse_bezier_patch.h"
#endif

#ifdef _MSC_VER
#define INLINE __forceinline
//define INLINE inline
#elif defined(__GNUC__)
#define INLINE __inline__
//define INLINE inline
#else
#define INLINE inline
#endif

namespace mallie{

#if USE_FLOAT 
#define REAL float
#define vector2x vector2f
#define vector3x vector3f
#define matrix3x matrix3f

static const REAL EPS = 1e-4;
static const REAL EPS2 = std::numeric_limits<float>::min();
static const REAL UVEPS = 1.0/32.0;
static const int  MAX_LEVEL = 10;
    
#else 
#define REAL double
#define vector2x vector2d
#define vector3x vector3d
#define matrix3x matrix3d

static const REAL EPS = 1e-5;
static const REAL EPS2 = std::numeric_limits<double>::min();
static const REAL UVEPS = 1.0/32.0;
static const int  MAX_LEVEL = 10;

#endif

    static
    inline vector3x Conv(const real3& r)
    {
        return vector3x(r[0], r[1], r[2]); 
    }
    
    static
    inline real3 Conv(const vector3x& r)
    {
        return real3(r[0], r[1], r[2]); 
    }

    
    struct range_AABB{
        float tmin;
        float tmax;
    };

    static inline bool intersectAABB(range_AABB* rng, const vector3 & min, const vector3 & max, const Ray & r, float tmin, float tmax)
    {
        //int phase = r.phase();
        int sign[3] = {r.dirSign[0],r.dirSign[1],r.dirSign[2]};
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
    
    template<class V, class PATCH>
    static void get_minmax(V& min, V& max, const PATCH& patch)
    {
    	static const REAL EPS3 = EPS*1e-3;
    	patch.minmax(min,max);
        //min = patch.min();
        //max = patch.max(); 
        min -= V(EPS3,EPS3,EPS3);
        max += V(EPS3,EPS3,EPS3);
    }
    
    template<class V>
    static
    void get_minmax(V& min, V& max, const std::vector<V>& p)
    {
        min = max = p[0];
        for(size_t i = 1;i<p.size();i++)
        {
            for(int j=0;j<3;j++)
            {
                min[j] = (min[j]>p[i][j])?p[i][j]:min[j];//if(min[j]>p[i][j])min[j]=p[i][j];
                max[j] = (max[j]<p[i][j])?p[i][j]:max[j];//if(max[j]<p[i][j])max[j]=p[i][j];
            }
        }
    }
    
    static 
    void getZAlign(matrix3& mat, const vector3x &dir)
    {
        vector3 z = dir;
        int plane = 0;
        if(fabs(z[1])<fabs(z[plane]))plane = 1;
        if(fabs(z[2])<fabs(z[plane]))plane = 2;
        vector3 x = vector3(0,0,0);
        x[plane] = 1;
        vector3 y = normalize(cross(z,x));
        x = cross(y,z);

        mat = matrix3(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2]);
    }
    
    static 
    void getZAlign(matrix4& mat, const Ray &r)
    {
        vector3x org = Conv(r.org);
        vector3x dir = Conv(r.dir);

        vector3x z = dir;
        int plane = 0;
        if(fabs(z[1])<fabs(z[plane]))plane = 1;
        if(fabs(z[2])<fabs(z[plane]))plane = 2;
        vector3x x = vector3(0,0,0);
        x[plane] = 1;
        vector3x y = normalize(cross(z,x));
        x = cross(y,z);
        matrix4 rot = matrix4(
                              x[0],x[1],x[2],0,
                              y[0],y[1],y[2],0,
                              z[0],z[1],z[2],0,
                              0,0,0,1
                              );
        matrix4 trs = matrix4(
                              1,0,0,-org[0],
                              0,1,0,-org[1],
                              0,0,1,-org[2],
                              0,0,0,1
                              );
        mat = rot*trs;
    }
    
    struct uvt_info
    {
        REAL u;
        REAL v;
        REAL t;
    };
    
    template<class VECTOR>
    static INLINE bool scan_minmax(const VECTOR& p)
    {
        size_t sz = p.size();
        int nUpper = 0;
        int nLower = 0;
        for(size_t i = 0;i<sz;i++){
            if(p[i][1]>0)     nUpper++;
            else nLower++;
            if(nUpper&&nLower)return true;
        }
        return false;
    }

    static 
    inline void bezier_minmax_(vector2x& min, vector2x& max, const vector2x p[], int n){
        min = max = p[0];
        for(int i = 1;i<n;i++){
            for(int j = 0;j<2;j++){
                min[j] = (min[j]>p[i][j])?p[i][j]:min[j];
                max[j] = (max[j]<p[i][j])?p[i][j]:max[j];
            }
        }
    }
    static 
    inline void bezier_minmax_(vector3x& min, vector3x& max, const vector3x p[], int n){
        min = max = p[0];
        for(int i = 1;i<n;i++){
            for(int j = 0;j<3;j++){
                min[j] = (min[j]>p[i][j])?p[i][j]:min[j];
                max[j] = (max[j]<p[i][j])?p[i][j]:max[j];
            }
        }
    }

    
    template<class T, int Sz>
    class static_vector
    {
    public:
        static_vector():sz_(0){}
        static_vector(int sz):sz_(sz){}
        size_t size()const{return sz_;}
        void push_back(const T& v){e_[sz_] = v;sz_++;}
        void pop_back(){sz_--;}
              T& operator[](size_t i)     {return e_[i];}
        const T& operator[](size_t i)const{return e_[i];}
    private:
        T e_[Sz];
        int sz_;
    };

	template<class T, size_t Sz=64>
	class static_bezier_patch
	{
	public:
		typedef T value_type;
		typedef static_bezier_patch<T,Sz> this_type;

		static const int DEFAULT_ORDER = 4;
	public:
		static_bezier_patch()
			:nu_(DEFAULT_ORDER), nv_(DEFAULT_ORDER){}
		static_bezier_patch(const bezier_patch<T>& rhs)
			:nu_(rhs.get_nu()), nv_(rhs.get_nv())
		{
			int nu = nu_;
			int nv = nv_;
			assert(nu*nv<=Sz);
			for(int j=0;j<nv;j++){
				for(int i=0;i<nu;i++){
					cp_[j*nu+i] = rhs.get_cp_at(i,j);
				}
			}
		}

		T evaluate     (float u, float v)const
		{
			return bezier_evaluate( &(cp_[0]), nu_, nv_, u, v);
		}

		int get_nu()const{return nu_;}
		int get_nv()const{return nv_;}
		T get_at(int i, int j)const{return cp_[nu_*j+i];}
		T get_cp_at(int i, int j)const{return get_at(i,j);}

/*
		T min()const
		{
			return bezier_min( &(cp_[0]), nu_*nv_ );
		}

		T max()const
		{
			return bezier_max( &(cp_[0]), nu_*nv_ );
		}
*/
        void minmax(T& min, T& max)const
        {
            bezier_minmax_(min,max, &(cp_[0]), nu_*nv_ );
        }
		this_type& transform(const matrix3x& m)
		{
			size_t sz = nu_*nv_;
			for(size_t i = 0;i<sz;i++){
				cp_[i] = m*cp_[i];
			}
			return *this;
		}

		this_type& transform(const matrix4& m)
		{
			size_t sz = nu_*nv_;
			for(size_t i = 0;i<sz;i++){
				cp_[i] = m*cp_[i];
			}
			return *this;
		}
		public:
		void split(this_type patches[4], float u, float v)const
		{
			int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			for(int i=0;i<4;i++){
				patches[i].nu_ = nu;
				patches[i].nv_ = nv;	
				//patches[i].cp_.resize(sz);
			}
			bezier_split(
				&(patches[0].cp_[0]),
				&(patches[1].cp_[0]),
				&(patches[2].cp_[0]),
				&(patches[3].cp_[0]),
				&(cp_[0]),
				nu,nv,
				u,v
			);
		}
        void split_u(this_type patches[2], float u)const
        {
            int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			for(int i=0;i<2;i++){
				patches[i].nu_ = nu;
				patches[i].nv_ = nv;	
				//patches[i].cp_.resize(sz);
			}
			bezier_split_u(
                &(patches[0].cp_[0]),
                &(patches[1].cp_[0]),
                &(cp_[0]),
                nu,nv,
                u
            );
        }
        void split_v(this_type patches[2], float v)const
        {
            int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			for(int i=0;i<2;i++){
				patches[i].nu_ = nu;
				patches[i].nv_ = nv;	
				//patches[i].cp_.resize(sz);
			}
			bezier_split_v(
                &(patches[0].cp_[0]),
                &(patches[1].cp_[0]),
                &(cp_[0]),
                nu,nv,
                v
            );
        }

		void crop_u(this_type& patch, float u0, float u1)const
        {
			int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			{
				patch.nu_ = nu;
				patch.nv_ = nv;	
				//patch.cp_.resize(sz);
			}
			bezier_crop_u(
				&(patch.cp_[0]),
                &(cp_[0]),
                nu,nv,	
				u0,u1
			);
		}

		void crop_v(this_type& patch, float v0, float v1)const
        {
			int nu = nu_;
			int nv = nv_;
			int sz = nu*nv;
			{
				patch.nu_ = nu;
				patch.nv_ = nv;	
				//patch.cp_.resize(sz);
			}
			bezier_crop_v(
				&(patch.cp_[0]),
                &(cp_[0]),
                nu,nv,	
				v0,v1
			);
		}
	private:
		int nu_;
		int nv_;
		T cp_[Sz];
	};
    
//define converx intersection 
#if 0  
    static
    INLINE float cross(float u1, float u2, float v1, float v2)
    {
        return -v1*u2+u1*v2;
    }
    
    static
    INLINE float check_line(const vector2x& a, const vector2x& b, const vector2x& p)
    {
        vector2x ab = b-a;
        vector2x ap = p-a;
        
        return cross(ab[0],ab[1],ap[0],ap[1]);
    }
    
    template<class VECTOR>
    static INLINE void push_upper_chain(VECTOR& chain, const vector2x& p)
    {
        int sz = (int)chain.size();
        if(sz < 2){
            chain.push_back(p);
        }else{
            vector2x a = chain[sz-2];//
            vector2x b = chain[sz-1];//
            if(0<=check_line(a,b,p)){
                chain.pop_back();
                push_upper_chain(chain,p);
            }else{
                chain.push_back(p);
            }
        }
    }
    
    template<class VECTOR>
    static INLINE void push_lower_chain(VECTOR& chain, const vector2x& p)
    {
        int sz = (int)chain.size();
        if(sz < 2){
            chain.push_back(p);
        }else{
            vector2x a = chain[sz-2];//
            vector2x b = chain[sz-1];//
            if(0>=check_line(a,b,p)){
                chain.pop_back();
                push_lower_chain(chain,p);
            }else{
                chain.push_back(p);
            }
        }
    }
        
    template<class VECTOR>
    static INLINE void convex_hull(VECTOR& upper_chain, VECTOR& lower_chain, const std::vector<vector2x>& points) 
    {
        using namespace std;
        int sz = (int)points.size();
        if(sz<2)return;
        upper_chain.push_back(points[0]);
        upper_chain.push_back(points[1]);
        for(int i = 2;i<sz;i++){
            push_upper_chain(upper_chain,points[i]);
        }
        lower_chain.push_back(points[0]);
        lower_chain.push_back(points[1]);
        for(int i = 2;i<sz;i++){
            push_lower_chain(lower_chain,points[i]);
        }
    }
    
    template<class VECTOR>
    static INLINE int scan_chain(float t[], const VECTOR& chain){
        int n = 0;
        size_t sz = chain.size();
        if(sz<1)return 0;
        for(size_t i = 0;i<sz-1;i++)
        {
            if(chain[i][1]*chain[i+1][1]<=0)
            {
                float a = fabs(chain[i  ][1]);
                float b = fabs(chain[i+1][1]);
                float h = a+b;
                if(h){
                    float k = chain[i][0] + a*(chain[i+1][0]-chain[i][0])/h;
                    t[n] = k;
                    n++;
                }
            }
        }
        return n;
    }
    
    template<class VECTOR>
    static int scan_convex(float ts[], const VECTOR& p)
    {
        if(!scan_minmax(p))
        {          
            return 0;
        }

        size_t sz = p.size();
        if(sz<=4){
            static_vector<vector2x,4> upper_chain;
            static_vector<vector2x,4> lower_chain;
            convex_hull(upper_chain,lower_chain,p);
            
            int n = 0;
            n += scan_chain(&ts[n],upper_chain);
            n += scan_chain(&ts[n],lower_chain);
            return n;
        }else{
            std::vector<vector2x> upper_chain;upper_chain.reserve(sz);
            std::vector<vector2x> lower_chain;lower_chain.reserve(sz);
            convex_hull(upper_chain,lower_chain,p);
            
            int n = 0;
            n += scan_chain(&ts[n],upper_chain);
            n += scan_chain(&ts[n],lower_chain);
            return n;
        }
    }
#else   
    static 
    INLINE float slope(const vector2x& a, const vector2x& b)
    {
    	vector2x dif = b-a;
    	return fabs(dif[0]/dif[1]);
    }
    
    static
   	INLINE float dist(const vector2x& p0, const vector2x& p1)
   	{
		float a = fabs(p0[1]);
		float b = fabs(p1[1]);
		float h = a+b;
		float t = p0[0] + a*(p1[0]-p0[0])/h;
		return t;
   	}
    
    template<class VECTOR>
    static INLINE int scan_convex(REAL ts[], const VECTOR& p)
    {
        if(!scan_minmax(p))
        {          
            return 0;
        }
        int n = 0;
        {
	        int k = -1;
	        int l = -1;
	        REAL current = std::numeric_limits<REAL>::max();
	        int sz = (int)p.size();
	        for(int i = 1;i<sz;i++)
	        {
	        	if(p[i][1]*p[0][1]<0){
	        		float s = slope(p[i], p[0]);
	        		if(s<current){
	        			current = s;
	        			k = i;
	        		}
	        	}
	        }
	        if(k<0)return 0;
	        current = 0;
	        for(int i = 0;i<k;i++)
	        {
	        	if(p[i][1]*p[k][1]<0){
	        		float s = slope(p[i], p[k]);
	        		if(current<s){
	        			current = s;
	        			l = i;
	        		}
	        	}
	        }
	        if(l<0)return 0;
            ts[n++] = dist(p[l],p[k]);
        }
        {
	        int k = -1;
	        int l = -1;
	        REAL current = std::numeric_limits<REAL>::max();
	        int sz = (int)p.size();
	        for(int i = 0;i<sz-1;i++)
	        {
	        	if(p[i][1]*p[sz-1][1]<0){
	        		float s = slope(p[i], p[sz-1]);
	        		if(s<current){
	        			current = s;
	        			k = i;
	        		}
	        	}
	        }
	        if(k<0)return 0;
	        current = 0;
	        for(int i = k+1;i<sz;i++)
	        {
	        	if(p[i][1]*p[k][1]<0){
	        		float s = slope(p[i], p[k]);
	        		if(current<s){
	        			current = s;
	        			l = i;
	        		}
	        	}
	        }
	        if(l<0)return 0;
            ts[n++] = dist(p[k],p[l]);
        }
        
     	return n;   
    }
#endif
    

    static
    inline REAL lerp(REAL a, REAL b, REAL t)
    {
        return a+(b-a)*t;
    }

	template<class PATCH>
    static void crop_u(PATCH& out, const PATCH& patch, REAL u0, REAL u1)
    {
#if 1
        patch.crop_u(out, u0, u1);
#else
        assert(u0<u1);
        PATCH a[2];
        PATCH b[2];
        patch.split_u(a,u1);
        a[0].split_u(b, u0/u1);
        out.swap(b[1]);
#endif
    }

	template<class PATCH>
    static void crop_v(PATCH& out, const PATCH& patch, REAL v0, REAL v1)
    {
#if 1
        patch.crop_v(out, v0, v1);
#else
        assert(v0<v1);
        PATCH a[2];
        PATCH b[2];
        patch.split_v(a,v1);
        a[0].split_v(b, v0/v1);
        out.swap(b[1]);
#endif
    }

	template<class VECTOR>
    static bool x_check(REAL rng[2], const VECTOR& curve)
    {
        REAL t[2] = {0};
        int nn = scan_convex(t, curve);
        if(nn){
            if(nn != 2)return false;

            float t0 = t[0];
            float t1 = t[1];

            if(t0>t1)std::swap(t0,t1);

            rng[0] = t0;
            rng[1] = t1;
            return true;
        }
        return false;
    }

//#define PRINT_A() std::cout << "+"
//#define PRINT_B() std::cout << "-"
#define PRINT_A()
#define PRINT_B()

    static 
    inline vector2x normalize2(vector3x& v)
    {
        float l = v[0]*v[0]+v[1]*v[1];
        l = 1.0/sqrt(l);
        v[0]*=l;
        v[1]*=l;
        return vector2x(v[0],v[1]);
    }
    static
    inline vector3x cross2(const vector3x& x)
    {
        return vector3x(
                -x[1], //xyzzy
                +x[0], //yzxxz
                0
            );
    }
    
    static
    void getRotate(matrix3x& mat, const vector3x &dx)
    {
        static const vector3x z = vector3x(0,0,1);
        vector3x x = dx;
        x[2] = 0;
        normalize2(x);//x = normalize(x)
        vector3x y = cross2(x);//cross(z,x)
        mat = matrix3x(x[0],x[1],x[2],y[0],y[1],y[2],z[0],z[1],z[2]);
    }

    template<class PATCH>
    static inline vector3x getLv(const PATCH& patch)
    {
        int nu = patch.get_nu();
        int nv = patch.get_nv();
        
        //int c = nu>>1;
        //return patch.get_cp_at(c,   nv-1)-patch.get_cp_at(c,   0);

        return
            patch.get_cp_at(0,   nv-1)-patch.get_cp_at(0,   0)+
            patch.get_cp_at(nu-1,nv-1)-patch.get_cp_at(nu-1,0);
    }

	template<class PATCH>
    static inline vector3x getLu(const PATCH& patch)
    {
        int nu = patch.get_nu();
        int nv = patch.get_nv();
        
        //int c = nv>>1;
        //return patch.get_cp_at(nu-1,   c)-patch.get_cp_at(0,   c);

        return
            patch.get_cp_at(nu-1,   0)-patch.get_cp_at(0,   0)+
            patch.get_cp_at(nu-1,nv-1)-patch.get_cp_at(0,nv-1);
    }

    static
    double tri_area(const vector2x& a, const vector2x& b, const vector2x& c)
    {
      return (a[0]-c[0])*(b[1]-c[1])  -  (a[1]-c[1])*(b[0]-c[0]);
    }

    /*
    //site:
    //Robust and Numerically Stable Bezier Clipping Method for Ray Tracing NURBS Surfaces
    //
	template<class PATCH>
    static
    vector3x getLv2(const PATCH& patch)
    {
        static const float K1 = 0.5;
        static const float K2 = sqrt(3.0)*0.5;
        static const vector2x O(0,0);

        vector2x lu = normalize2(getLu(patch));
        vector2x lv = normalize2(getLv(patch));

        float c = fabs(dot(lu, lv));

        if(c>0.5){
            vector2x M = 0.5*(lu+lv);
            vector2x N1 = vector2x(M[1],-M[0]);
            vector2x N2 = vector2x(-M[1],M[0]);

            if(tri_area(lu,O,M)<0){
                std::swap(N1,N2);
            }

            vector2x R = K1*N1 + K2*M;
            return vector3x(R[0],R[1],0);          
        }else{
            return vector3x(lv[0],lv[1],0);                  
        }      
    }

	template<class PATCH>
    static
    vector3x getLu2(const PATCH& patch)
    {
        static const float K1 = 0.5;
        static const float K2 = sqrt(3.0)*0.5;
        static const vector2x O(0,0);

        vector2x lu = normalize2(getLu(patch));
        vector2x lv = normalize2(getLv(patch));

        float c = fabs(dot(lu, lv));

        if(c>0.5){
            vector2x M = 0.5*(lu+lv);
            vector2x N1 = vector2x(M[1],-M[0]);
            vector2x N2 = vector2x(-M[1],M[0]);

            if(tri_area(lu,O,M)<0){
                std::swap(N1,N2);
            }

            vector2x R = K1*N2 + K2*M;
            return vector3x(R[0],R[1],0);          
        }else{
            return vector3x(lu[0],lu[1],0);            
        }      
    }
    */
	template<class PATCH>
    static void rotate_u(PATCH& patch)
    {
        vector3x dx = getLv(patch);
        matrix3x rot;
        getRotate(rot, dx);
        patch.transform(rot);
    }

	template<class PATCH>
    static void rotate_v(PATCH& patch)
    {
        vector3x dx = getLu(patch);
        matrix3x rot;
        getRotate(rot, dx);
        patch.transform(rot);
    }

	template<class PATCH>
    static bool get_range_u(REAL out[2], const PATCH& patch)
    {
        int nu = patch.get_nu();
        int nv = patch.get_nv();

        REAL t0 = 1;
        REAL t1 = 0;
        
        if(nu<=16){
	        static_vector<vector2x,16> curve(nu);
	        REAL delta = 1.0/(nu-1);
	        for(int i=0;i<nu;i++){
	            curve[i][0] = i*delta;
	        }
	        REAL rng[2];
	        for(int j=0;j<nv;j++){
	            for(int i=0;i<nu;i++){
	                curve[i][1] = (patch.get_cp_at(i,j))[1];
	            }
	
	            if(x_check(rng,curve)){
	                if(rng[1]-rng[0]<delta){
	                    t0 = std::min(t0,rng[0]);
	                    t1 = std::max(t1,rng[1]);
	                }else{
	                    return false;
	                }
	            }else{
	                return false;
	            }
	        }
        }else{
        	std::vector<vector2x> curve(nu);
	        REAL delta = 1.0/(nu-1);
	        for(int i=0;i<nu;i++){
	            curve[i][0] = i*delta;
	        }
	        REAL rng[2];
	        for(int j=0;j<nv;j++){
	            for(int i=0;i<nu;i++){
	                curve[i][1] = (patch.get_cp_at(i,j))[1];
	            }
	
	            if(x_check(rng,curve)){
	                if(rng[1]-rng[0]<delta){
	                    t0 = std::min(t0,rng[0]);
	                    t1 = std::max(t1,rng[1]);
	                }else{
	                    return false;
	                }
	            }else{
	                return false;
	            }
	        }
        }

        out[0] = t0;
        out[1] = t1;
        return true;
    }

	template<class PATCH>
    static bool get_range_v(REAL out[2], const PATCH& patch)
    {
        int nu = patch.get_nu();
        int nv = patch.get_nv();

        REAL t0 = 1;
        REAL t1 = 0;
        
        if(nv<=16){
	        static_vector<vector2x,16> curve(nv);
	        REAL delta = 1.0/(nv-1);
	        for(int i=0;i<nv;i++){
	            curve[i][0] = i*delta;
	        }
	
	        REAL rng[2];
	        for(int i=0;i<nu;i++){
	            for(int j=0;j<nv;j++){          
	                curve[j][1] = (patch.get_cp_at(i,j))[1];
	            }
	            if(x_check(rng,curve)){
	                if(rng[1]-rng[0]<delta){
	                    t0 = std::min(t0,rng[0]);
	                    t1 = std::max(t1,rng[1]);
	                }else{
	                    return false;
	                }
	            }else{
	                return false;
	            }
	        }
        }else{
        	std::vector<vector2x> curve(nv);
	        REAL delta = 1.0/(nv-1);
	        for(int i=0;i<nv;i++){
	            curve[i][0] = i*delta;
	        }
	
	        REAL rng[2];
	        for(int i=0;i<nu;i++){
	            for(int j=0;j<nv;j++){          
	                curve[j][1] = (patch.get_cp_at(i,j))[1];
	            }
	            if(x_check(rng,curve)){
	                if(rng[1]-rng[0]<delta){
	                    t0 = std::min(t0,rng[0]);
	                    t1 = std::max(t1,rng[1]);
	                }else{
	                    return false;
	                }
	            }else{
	                return false;
	            }
	        }
        	
        }

        out[0] = t0;
        out[1] = t1;
        return true;
    }

	template<class PATCH>
    static bool test_bezier_clip_l(const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax);
    
	template<class PATCH>
	static bool test_bezier_clip_u(const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps);
    template<class PATCH>
	static bool test_bezier_clip_v(const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps);

	template<class PATCH>
    static bool test_bezier_clip_l(uvt_info* info, const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax);
    template<class PATCH>
	static bool test_bezier_clip_u(uvt_info* info, const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps);
    template<class PATCH>
	static bool test_bezier_clip_v(uvt_info* info, const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps);

	template<class PATCH>
    static bool test_bezier_clip_l(const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax)
    {
        REAL t;
        REAL u;
        REAL v;
#if DIRECT_BILINEAR
        vector3x P[4];
        int nu = patch.get_nu();
        int nv = patch.get_nv();
        P[0] = patch.get_cp_at(0,0);
        P[1] = patch.get_cp_at(nu-1,0);
        P[2] = patch.get_cp_at(0,nv-1);
        P[3] = patch.get_cp_at(nu-1,nv-1);

        if(test_bilinear_patch(&t, &u, &v, P, zmin, zmax)){
            return true;
        }
        return false;
#else
        //bool bRet = false;
        int nPu = patch.get_nu()-1;
        int nPv = patch.get_nv()-1;

        REAL du = REAL(1)/nPu;
        REAL dv = REAL(1)/nPv;

        vector3x P[4];
        for(int j=0;j<nPv;j++)
        {
            for(int i=0;i<nPu;i++)
            {
                P[0]=patch.get_cp_at(i  ,j  );
                P[1]=patch.get_cp_at(i+1,j  );
                P[2]=patch.get_cp_at(i  ,j+1);
                P[3]=patch.get_cp_at(i+1,j+1);
                if(test_bilinear_patch(&t, &u, &v, P, zmin, zmax)){
                    return true;
                }
            }
        }
        return false;   
#endif
    }

	template<class PATCH>
    static bool test_bezier_clip_l(uvt_info* info, const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax)
    {
        REAL t;
        REAL u;
        REAL v;
#if DIRECT_BILINEAR
        vector3x P[4];
        int nu = patch.get_nu();
        int nv = patch.get_nv();
        P[0] = patch.get_cp_at(0,0);
        P[1] = patch.get_cp_at(nu-1,0);
        P[2] = patch.get_cp_at(0,nv-1);
        P[3] = patch.get_cp_at(nu-1,nv-1);

        if(test_bilinear_patch(&t, &u, &v, P, zmin, zmax)){
            u = u0*(1-u)+u1*u;
            v = v0*(1-v)+v1*v;
            info->u = u;
            info->v = v;
            info->t = t;
            return true;
        }
        return false;
#else
        bool bRet = false;
        int nPu = patch.get_nu()-1;
        int nPv = patch.get_nv()-1;

        REAL du = REAL(1)/nPu;
        REAL dv = REAL(1)/nPv;

        REAL uu = 0;
        REAL vv = 0;

        vector3x P[4];
        for(int j=0;j<nPv;j++)
        {
            for(int i=0;i<nPu;i++)
            {
                P[0]=patch.get_cp_at(i  ,j  );
                P[1]=patch.get_cp_at(i+1,j  );
                P[2]=patch.get_cp_at(i  ,j+1);
                P[3]=patch.get_cp_at(i+1,j+1);
                if(test_bilinear_patch(&t, &u, &v, P, zmin, zmax)){

                    u = lerp(i*du,(i+1)*du,u);
                    v = lerp(j*dv,(j+1)*dv,v);

                    uu = u;
                    vv = v;

                    zmax = t;

                    u = lerp(u0,u1,u);
                    v = lerp(v0,v1,v);

                    info->u = u;
                    info->v = v;
                    info->t = t;
                    bRet = true;
                }
            }
        }

        if(bRet)
        {
            vector3x p = patch.evaluate(uu,vv);
            info->t = p[2];
        }

        return bRet;    
#endif
    }

    static 
    inline bool is_eps(const vector3x& min, const vector3x& max, float eps)
    {
        //float zw = max[2]-min[2];
        //if(zw<=eps)return true;
        
        REAL xd = std::max<REAL>(fabs(min[0]),max[0]);
        REAL yd = std::max<REAL>(fabs(min[1]),max[1]);
        if(!(xd<=eps&&yd<=eps))return false;

        REAL xw = max[0]-min[0];
        REAL yw = max[1]-min[1];
        //
        if(xw<=eps && yw<=eps)return true;
        else return false;
    }

    static
    inline bool is_level(int level, int max_level)
    {
        return (level>=max_level);
        //return false;
    }

    static
    inline REAL log2(REAL x)
    {
        static const REAL LOG2 = 1.0/log(2.0);
        return log(x)*LOG2;    
    }

    static
    inline bool is_clip(int level, int nc)
    {
        static const int div_level[] = 
        {
            1,
            1,1,
            2,2,//4
            3,3,3,3,//8
            4,4,4,4,4,4,4,4,//16
        };
        int nlevel = 0; 
        if(nc<=16){
            nlevel = div_level[nc];
        }else{
            nlevel = (int)ceil(log2(nc));
        }
        return nlevel*2<=level;
        //return true;
    }

	template<class PATCH>
    static void coarse_sort(int order[2], PATCH tmp[2])
    {
#if USE_COARSESORT
        REAL zs[2];
        for(int i=0;i<2;i++){
            int u = tmp[i].get_nu()>>1;
            int v = tmp[i].get_nv()>>1;

            zs[i] = tmp[i].get_cp_at(u,v)[2];
        }
        if(zs[0]<=zs[1]){   
            order[0] = 0;
            order[1] = 1;
        }else{
            order[0] = 1;
            order[1] = 0;
        }
#endif
    }

	template<class PATCH>
    static bool test_bezier_clip_u( const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps)
    {
        PATCH tpatch(patch);
        rotate_u(tpatch);

        vector3x min;
        vector3x max;
        get_minmax(min, max, tpatch);
        if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<zmin || zmax<min[2])return false;//z
        //if(max[0]-min[0]<=EPS2)return false;
        if(max[1]-min[1]<=EPS2)return false;

        if( is_eps(min,max,eps) || is_level(level,max_level) ){
            return test_bezier_clip_l(patch, u0, u1, v0, v1, zmin, zmax);
        }else{
            REAL tw = 1;
            REAL tt0 = 1;
            REAL tt1 = 0;
#if USE_BEZIERCLIP
            if(is_clip(level,tpatch.get_nu())){
                REAL rng[2];
                if(get_range_u(rng, tpatch)){
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }
#endif
            if(tw>=0.4){
                PATCH tmp[2];
                patch.split_u(tmp, 0.5);
                REAL um = (u0+u1)*0.5;
                REAL ut[] = {u0,um,um,u1};

                int order[2]={0,1};
                coarse_sort(order, tmp);
                //bool bRet = false;
                for(int i=0;i<2;i++){
                    if(test_bezier_clip_v(tmp[order[i]], ut[2*order[i]], ut[2*order[i]+1], v0, v1, zmin, zmax, level+1, max_level, eps)){
                        return true;
                    }
                }
                return false;
            }else{
                tt0 = std::max<REAL>(0.0,tt0-UVEPS);
                tt1 = std::min<REAL>(tt1+UVEPS,1.0);
                REAL ut[] = {lerp(u0,u1,tt0),lerp(u0,u1,tt1)};
                PATCH tmp;
                crop_u(tmp, patch, tt0, tt1);
                return test_bezier_clip_v(tmp, ut[0], ut[1], v0, v1, zmin, zmax, level+1, max_level, eps);
            }
        }
        return false;
    }

	template<class PATCH>
    static bool test_bezier_clip_v(const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps)
    {
        PATCH tpatch(patch);
        rotate_v(tpatch);

        vector3x min;
        vector3x max;
        get_minmax(min, max, tpatch);
        if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<zmin || zmax<min[2])return false;//z
        //if(max[0]-min[0]<=EPS2)return false;
        if(max[1]-min[1]<=EPS2)return false;

        if( is_eps(min,max,eps) || is_level(level,max_level) ){
            return test_bezier_clip_l(patch, u0, u1, v0, v1, zmin, zmax);
        }else{
            REAL tw = 1;
            REAL tt0 = 1;
            REAL tt1 = 0;
#if USE_BEZIERCLIP
            if(is_clip(level,tpatch.get_nv())){
                REAL rng[2];
                if(get_range_v(rng, tpatch)){
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }
#endif
            
            if(tw>=0.4){
                PATCH tmp[2];
                patch.split_v(tmp, 0.5);
                REAL vm = (v0+v1)*0.5;
                REAL vt[] = {v0,vm,vm,v1};

                int order[2]={0,1};
                coarse_sort(order, tmp);
                //bool bRet = false;
                for(int i=0;i<2;i++){
                    if(test_bezier_clip_u(tmp[order[i]], u0, u1, vt[2*order[i]], vt[2*order[i]+1], zmin, zmax, level+1, max_level, eps)){
                        return true;
                    }
                }
                return false;
            }else{
                tt0 = std::max<REAL>(0.0,tt0-UVEPS);
                tt1 = std::min<REAL>(tt1+UVEPS,1.0);
                REAL vt[] = {lerp(v0,v1,tt0),lerp(v0,v1,tt1)};
                PATCH tmp;
                crop_v(tmp, patch, tt0, tt1);
                return test_bezier_clip_u(tmp, u0, u1, vt[0], vt[1], zmin, zmax, level+1, max_level, eps);
            }
        }
        return false;
    }

	template<class PATCH>
    static bool test_bezier_clip_u(uvt_info* info, const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps)
    {
        PATCH tpatch(patch);
        rotate_u(tpatch);

        vector3x min;
        vector3x max;
        get_minmax(min, max, tpatch);
        if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<zmin || zmax<min[2])return false;//z
        //if(max[0]-min[0]<=EPS2)return false;
        if(max[1]-min[1]<=EPS2)return false;

        if( is_eps(min,max,eps) || is_level(level,max_level) ){
            return test_bezier_clip_l(info, patch, u0, u1, v0, v1, zmin, zmax);
        }else{
            REAL tw = 1;
            REAL tt0 = 1;
            REAL tt1 = 0;
#if USE_BEZIERCLIP
            if(is_clip(level,tpatch.get_nu())){
                REAL rng[2];
                if(get_range_u(rng, tpatch)){
                    tt0 = rng[0];
                    tt1 = rng[1];                  
                    tw = tt1-tt0;
                }
            }
#endif
            
            if(tw>=0.4){
                PATCH tmp[2];
                patch.split_u(tmp, 0.5);
                REAL um = (u0+u1)*0.5;
                REAL ut[] = {u0,um,um,u1};

                int order[2]={0,1};
                coarse_sort(order, tmp);
                bool bRet = false;
                for(int i=0;i<2;i++){
                    if(test_bezier_clip_v(info, tmp[order[i]], ut[2*order[i]], ut[2*order[i]+1], v0, v1, zmin, zmax, level+1, max_level, eps)){
                        zmax = info->t;
                        bRet = true;
                    }
                }
                return bRet;
            }else{
                tt0 = std::max<REAL>(0.0,tt0-UVEPS);
                tt1 = std::min<REAL>(tt1+UVEPS,1.0);
                REAL ut[] = {lerp(u0,u1,tt0),lerp(u0,u1,tt1)};
                PATCH tmp;
                crop_u(tmp, patch, tt0, tt1);
                return test_bezier_clip_v(info, tmp, ut[0], ut[1], v0, v1, zmin, zmax, level+1, max_level, eps);
            }
        }
        return false;
    }

	template<class PATCH>
    static bool test_bezier_clip_v(uvt_info* info, const PATCH& patch, REAL u0, REAL u1, REAL v0, REAL v1, REAL zmin, REAL zmax, int level, int max_level, REAL eps)
    {
        PATCH tpatch(patch);
        rotate_v(tpatch);

        vector3x min;
        vector3x max;
        get_minmax(min, max, tpatch);
        if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<zmin || zmax<min[2])return false;//z
        //if(max[0]-min[0]<=EPS2)return false;
        if(max[1]-min[1]<=EPS2)return false;

        if( is_eps(min,max,eps) || is_level(level,max_level) ){
            return test_bezier_clip_l(info, patch, u0, u1, v0, v1, zmin, zmax);
        }else{
            REAL tw = 1;
            REAL tt0 = 1;
            REAL tt1 = 0;
#if USE_BEZIERCLIP
            if(is_clip(level,tpatch.get_nv())){
                REAL rng[2];
                if(get_range_v(rng, tpatch)){
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }
#endif
            
            if(tw>=0.4){
                PATCH tmp[2];
                patch.split_v(tmp, 0.5);
                REAL vm = (v0+v1)*0.5;
                REAL vt[] = {v0,vm,vm,v1};

                int order[2]={0,1};
                coarse_sort(order, tmp);
                bool bRet = false;
                for(int i=0;i<2;i++){
                    if(test_bezier_clip_u(info, tmp[order[i]], u0, u1, vt[2*order[i]], vt[2*order[i]+1], zmin, zmax, level+1, max_level, eps)){
                        zmax = info->t;
                        bRet = true;
                    }
                }
                return bRet;
            }else{
                tt0 = std::max<REAL>(0.0,tt0-UVEPS);
                tt1 = std::min<REAL>(tt1+UVEPS,1.0);
                REAL vt[] = {lerp(v0,v1,tt0),lerp(v0,v1,tt1)};
                PATCH tmp;
                crop_v(tmp, patch, tt0, tt1);
                return test_bezier_clip_u(info, tmp, u0, u1, vt[0], vt[1], zmin, zmax, level+1, max_level, eps);
            }
        }
        return false;
    }

	template<class PATCH>
    static bool test_bezier_patch(const PATCH& patch, REAL zmin, REAL zmax, REAL eps)
    {
        vector3 min;
        vector3 max;
        get_minmax(min, max, patch);
        if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<zmin || zmax<min[2])return false;//z

        return test_bezier_clip_u(patch, 0, 1, 0, 1, zmin, zmax, 0, MAX_LEVEL, eps);
    }
    
	template<class PATCH>
    static bool test_bezier_patch(uvt_info* info, const PATCH& patch, REAL zmin, REAL zmax, REAL eps)
    {
        vector3 min;
        vector3 max;
        get_minmax(min, max, patch);
        if(0<min[0] || max[0]<0)return false;//x
        if(0<min[1] || max[1]<0)return false;//y
        if(max[2]<zmin || zmax<min[2])return false;//z

        return test_bezier_clip_u(info, patch, 0, 1, 0, 1, zmin, zmax, 0, MAX_LEVEL, eps);
    }
    
    //-----------------------------------------------------------------------------
    
    static
    float get_eps(const vector3& min, const vector3& max)
    {
        //vector3 w = max-min;
        //float ww = std::max(w[0],std::max(w[1],w[2]));
        return EPS;//std::max(ww/4000.0,EPS);
    }

    bezier_patch_intersection::bezier_patch_intersection(const bezier_patch<vector3>& patch)
        :patch_(patch)
    {
        urange_[0] = 0;
        urange_[1] = 1;

        vrange_[0] = 0;
        vrange_[1] = 1;
        get_minmax(min_,max_, patch);
        min_ += -vector3(2*EPS,2*EPS,2*EPS);
        max_ += +vector3(2*EPS,2*EPS,2*EPS);

        eps_ = get_eps(min_,max_);
    }
    bezier_patch_intersection::bezier_patch_intersection(const bezier_patch<vector3>& patch, float u0, float u1, float v0, float v1)
        :patch_(patch)
    {
        float EPSILON = EPS*1024;
        urange_[0] = u0;
        urange_[1] = u1;
        vrange_[0] = v0;
        vrange_[1] = v1;
        get_minmax(min_,max_, patch);
        min_ += -vector3(2*EPS,2*EPS,2*EPS);
        max_ += +vector3(2*EPS,2*EPS,2*EPS);
        
        eps_ = get_eps(min_,max_);
    }
    
    bezier_patch_intersection::~bezier_patch_intersection()
    {
        ;
    }

    
    
    bool bezier_patch_intersection::test_internal    (                 const Ray& r, float tmin, float tmax)const
    {
        matrix4 mat;
        getZAlign(mat, r);

		int nu = patch_.get_nu();
		int nv = patch_.get_nv();

		if(nu*nv<=16){
            static_bezier_patch<vector3x,16> patch(patch_);
			patch.transform(mat);
			return test_bezier_patch(patch, tmin, tmax, eps_);
		}else{
			bezier_patch<vector3x> patch(patch_);
			patch.transform(mat);
			return test_bezier_patch(patch, tmin, tmax, eps_);
		}
    }
    
    bool bezier_patch_intersection::test_internal    (Intersection* info, const Ray& r, float tmin, float tmax)const
    {
        matrix4 mat;
        getZAlign(mat, r);

		int nu = patch_.get_nu();
		int nv = patch_.get_nv();

		bool bRet = false;
		uvt_info uvt;
		if(nu*nv<=16){
            static_bezier_patch<vector3x,16> patch(patch_);
			patch.transform(mat);
			bRet = test_bezier_patch(&uvt, patch, tmin, tmax, eps_);
		}else{
			bezier_patch<vector3x> patch(patch_);
			patch.transform(mat);
			bRet = test_bezier_patch(&uvt, patch, tmin, tmax, eps_);
		}

        if(bRet)
        {
            
            float t = uvt.t;
            float u = uvt.u;
            float v = uvt.v;      

            u = urange_[0]*(1-u) + urange_[1]*u;//global
            v = vrange_[0]*(1-v) + vrange_[1]*v;//global
            info->t = t;
            info->position = Conv(Conv(r.org) + t*Conv(r.dir));
            info->u = u;
            info->v = v;
            {
				//u = std::max<float>(0.01, std::min<float>(u,0.99));
				//v = std::max<float>(0.01, std::min<float>(v,0.99));
		
		        vector3 U = normalize(patch_.evaluate_deriv_u(u,v));
		        vector3 V = normalize(patch_.evaluate_deriv_v(u,v));
		        vector3 N = normalize(cross(U,V));
		        info->normal = Conv(N);
		        info->geometricNormal = Conv(N);
		        info->tangent  = Conv(U);
		        info->binormal = Conv(V);
            }
            //uvt_info *fsp = reinterpret_cast<uvt_info*>(info->freearea);
            //*fsp = uvt;
			//this->finalize(info, r, t);
            //info->p_intersection = this;
            return true;
        }
        return false;
    }
    
    bool bezier_patch_intersection::test    (                 const Ray& r, float tmin, float tmax)const
    {
        range_AABB rng;
        if(intersectAABB(&rng, min_, max_, r, tmin, tmax)){
            tmin = std::max<float>(tmin, rng.tmin);
            tmax = std::min<float>(tmax, rng.tmax);
            return this->test_internal(     r, tmin, tmax);
        }
        return false;
    }
    
    bool bezier_patch_intersection::test    (Intersection* info, const Ray& r, float tmin, float tmax)const
    {
        range_AABB rng;
        if(intersectAABB(&rng, min_, max_, r, tmin, tmax)){

            tmin = std::max<float>(tmin, rng.tmin);
            tmax = std::min<float>(tmax, rng.tmax);
            return this->test_internal(info, r, tmin, tmax);
        }
        return false;
    }
    
    bool bezier_patch_intersection::test    (                 const Ray& r, float dist)const
    {
        return this->test(      r, 0, dist);
    }
    
    bool bezier_patch_intersection::test    (Intersection* info, const Ray& r, float dist)const
    {
        return this->test(info, r, 0, dist);
    }
    
    void bezier_patch_intersection::finalize(Intersection* info, const Ray& r, float dist)const
    {
        info->position = Conv(Conv(r.org) + dist*Conv(r.dir));
        /*
        uvt_info uvt = *reinterpret_cast<uvt_info*>(info->freearea);
        
        float u = uvt.u;
        float v = uvt.v;

		//u = std::max<float>(0.01, std::min<float>(u,0.99));
		//v = std::max<float>(0.01, std::min<float>(v,0.99));

        vector3 U = normalize(patch_.evaluate_deriv_u(u,v));
        vector3 V = normalize(patch_.evaluate_deriv_v(u,v));
        vector3 N = cross(U,V);
        info->normal = N;
        info->geometric = N;
        info->tangent  = U;
        info->binormal = V;
        */
    }
    
    inline vector3 bezier_patch_intersection::min()const
    {
        return min_;
    }
    
    inline vector3 bezier_patch_intersection::max()const
    {
        return max_;
    }
        
}

