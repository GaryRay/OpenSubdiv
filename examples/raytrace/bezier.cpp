#include "bezier.h"

#include <vector>

#ifndef _WIN32
#include <alloca.h>
#endif

static const int BINOMIAL_TABLE[12][12] = {
    {1},                                        /*0*/
    {1,1},                                      /*1*/
    {1,2,1},                                    /*2*/ 
    {1,3,3,1},                                  /*3*/
    {1,4,6,4,1},                                /*4*/
    {1,5,10,10,5,1},                            /*5*/
    {1,6,15,20,15,6,1},                         /*6*/
    {1,7,21,35,35,21,7,1},                      /*7*/
    {1,8,28,56,70,56,28,8,1},                   /*8*/
    {1,9,36,84,126,126,84,36,9,1},              /*9*/
    {1,10,45,120,210,252,210,120,45,10,1},      /*10*/
    {1,11,55,165,330,462,462,330,165,55,11,1},  /*11*/
};

static
inline int binomial(int i, int j)
{
  assert(0<=i);
  assert(0<=j && j<=i);
  if(i<=11){
    return BINOMIAL_TABLE[i][j];
  }else{
    if(j==0||j==i){
      return 1;
    }else{
      return binomial(i-1,j-1)+binomial(i-1,j);
    }
  }
}

static
inline int sgn(int i, int j)
{
  return ((i+j)&1)?1:-1;
}

namespace {

	template<class T>
	class alloca_vector
	{
	public:
		alloca_vector(size_t sz)
		{
			ptr_ = (T*)alloca(sz*sizeof(T));
		}
		~alloca_vector()
		{
			;//
		}
		      T& operator[](size_t i){return ptr_[i];}
		const T& operator[](size_t i)const{return ptr_[i];}
	private:
		T* ptr_;
	};

}

namespace mallie{

	template<class FLOAT>
    static inline void bernstein_2(FLOAT t, FLOAT e[4]){
        FLOAT s = 1-t;
        
        e[0] = s;
        e[1] = t;
    }
    
	template<class FLOAT>
    static inline void bernstein_4(FLOAT t, FLOAT e[4]){
        FLOAT s = 1-t;
        FLOAT s2 = s*s;
        FLOAT t2 = t*t;

        e[0] = 1*s*s2;
        e[1] = 3*s2*t;
        e[2] = 3*s*t2;
        e[3] = 1*t*t2;
    }

	template<class FLOAT>
    static inline void bernstein_n(FLOAT t, FLOAT e[], int n){
		switch(n){
		case 2:
			bernstein_2(t,e);
			break;
		case 4:
			bernstein_4(t,e);
			break;
		default:
			{
				memset(e,0,sizeof(FLOAT)*n);//ZERO 
				FLOAT s = 1;
				for(int i = n-1;0<=i;i--){	//2,1,0
					int k = n-i;			//1,2,3
					for(int j = 0;j<k;j++){
						e[j] += sgn(i+n,j)*binomial(n-1-j,i)*s;
					}
					s *= t;
				}
				for(int i = 0;i<n;i++){
					e[i] *= (FLOAT)(binomial(n-1,i));
				}
			}
			break;
		}
    }

#define DEC_BEZIER_BERNSTEIN(FLOAT)                  \
    void bernstein(FLOAT t, FLOAT e[4]){             \
        bernstein_4(t, e);                           \
    }                                                \
    void bernstein(FLOAT t, FLOAT e[], int n){       \
        bernstein_n(t, e, n);                        \
    }    

	DEC_BEZIER_BERNSTEIN(float)
	DEC_BEZIER_BERNSTEIN(double)

#undef DEC_BEZIER_BERNSTEIN



    //-----------------------------------------------

    template<class FLOAT>
    static inline void bernstein_deriv_2(FLOAT t, FLOAT e[4]){
        e[0] = -1;
        e[1] = +1;
    }

	template<class FLOAT>
    static inline void bernstein_deriv_4(FLOAT t, FLOAT e[4]){
        FLOAT t2 = t*t;

        e[0] = 3*t2*-1 + 2*t* 3 + -3;
        e[1] = 3*t2* 3 + 2*t*-6 +  3;
        e[2] = 3*t2*-3 + 2*t* 3;
        e[3] = 3*t2* 1;
    }

	template<class FLOAT>
	static inline void bernstein_deriv_n(FLOAT t, FLOAT e[], int n){
		switch(n){
		case 2:
			bernstein_deriv_2(t,e);
			break;
		case 4:
			bernstein_deriv_4(t,e);
			break;
		default:
			{
				memset(e,0,sizeof(FLOAT)*n);//ZERO 
				FLOAT s = 1;
				for(int i = n-2;0<=i;i--){//
					int k = n-i;//
					for(int j = 0;j<k;j++){
						e[j] += sgn(i+n,j)*binomial(n-1-j,i)*s*(n-1-i);
					}
					s *= t;
				}
				for(int i = 0;i<n;i++){
					e[i] *= (FLOAT)(binomial(n-1,i));
				}
			}
			break;
		}
	}

#define DEC_BEZIER_BERNSTEIN_DERIV(FLOAT)            \
    void bernstein_deriv(FLOAT t, FLOAT e[4]){       \
        bernstein_deriv_4(t, e);                     \
    }                                                \
    void bernstein_deriv(FLOAT t, FLOAT e[], int n){ \
        bernstein_deriv_n(t, e, n);                  \
    }    

	DEC_BEZIER_BERNSTEIN_DERIV(float)
	DEC_BEZIER_BERNSTEIN_DERIV(double)

#undef DEC_BEZIER_BERNSTEIN_DERIV

    //-----------------------------------------------

#define DEC_BEZIER_MUL_COEF(T, FLOAT)                                                           \
    static inline T mul_coef(const FLOAT e[4], const T& a, const T& b, const T& c, const T& d){ \
        return e[0]*a + e[1]*b + e[2]*c + e[3]*d;                                               \
    }                                                                                           \
    static inline T mul_coef(const FLOAT e[], const T v[], int n){                              \
        T tmp = e[0]*v[0];                                                                      \
        for(int i = 1;i<n;i++)tmp += e[i]*v[i];                                                 \
        return tmp;                                                                             \
    }            

	DEC_BEZIER_MUL_COEF(float, float)
	DEC_BEZIER_MUL_COEF(double, double)

	//DEC_BEZIER_MUL_COEF(vector2f, float)
	DEC_BEZIER_MUL_COEF(vector3f, float)
	//DEC_BEZIER_MUL_COEF(vector4f, float)
	//DEC_BEZIER_MUL_COEF(vector2d, double)
	DEC_BEZIER_MUL_COEF(vector3d, double)
	//DEC_BEZIER_MUL_COEF(vector4d, double)

#undef DEC_BEZIER_MUL_COEF
  //-----------------------------------------------
	template<class T,class FLOAT>
    static inline T bezier_evaluate_2(const T p[2], FLOAT t){
        return p[0]+(p[1]-p[0])*t;
    }

    template<class T,class FLOAT>
    static inline T bezier_evaluate_4(const T p[4], FLOAT t){
        FLOAT B[4];
        bernstein_4(t,B);
        return mul_coef(B,p[0],p[1],p[2],p[3]);
    }

    template<class T,class FLOAT>
    static inline T bezier_evaluate_4(const T p[4], const FLOAT coe[4]){
        return mul_coef(coe,p[0],p[1],p[2],p[3]);
    }

    template<class T,class FLOAT>
    static inline T bezier_evaluate_n(const T p[], int n, FLOAT t){
		switch(n){
		case 2:
			return bezier_evaluate_2(p, t);
		case 4:
			return bezier_evaluate_4(p, t);
		default:
			{
				alloca_vector<FLOAT> B(n);
				bernstein_n(t,&B[0],n);
				return mul_coef(&B[0],p,n);
			}
		}
    }

#define DEC_BEZIER_LINE_EVALUATE(TYPE,FLOAT)                                                  \
  TYPE bezier_evaluate(const TYPE p[4], FLOAT t){return bezier_evaluate_4(p,t);}              \
  TYPE bezier_evaluate(const TYPE p[4], const FLOAT coe[4]){return bezier_evaluate_4(p,coe);} \
  TYPE bezier_evaluate(const TYPE p[] , int n, FLOAT t){return bezier_evaluate_n(p,n,t);}

    DEC_BEZIER_LINE_EVALUATE(float,float)
    //DEC_BEZIER_LINE_EVALUATE(vector2f,float)
    DEC_BEZIER_LINE_EVALUATE(vector3f,float)
    //DEC_BEZIER_LINE_EVALUATE(vector4f,float)
	DEC_BEZIER_LINE_EVALUATE(double,double)
    //DEC_BEZIER_LINE_EVALUATE(vector2d,double)
    DEC_BEZIER_LINE_EVALUATE(vector3d,double)
    //DEC_BEZIER_LINE_EVALUATE(vector4d,double)

#undef DEC_BEZIER_LINE_EVALUATE

	template<class T,class FLOAT>
    static inline T bezier_evaluate_deriv_2(const T p[2], FLOAT t){
		return p[1]-p[0];
    }

    template<class T,class FLOAT>
    static inline T bezier_evaluate_deriv_4(const T p[4], FLOAT t){
        FLOAT B[4];
        bernstein_deriv(t,B);
        return mul_coef(B,p[0],p[1],p[2],p[3]);
    }

    template<class T,class FLOAT>
    static inline T bezier_evaluate_deriv_4(const T p[4], const FLOAT coe[4]){
        return mul_coef(coe,p[0],p[1],p[2],p[3]);
    }

    template<class T,class FLOAT>
    static inline T bezier_evaluate_deriv_n(const T p[], int n, FLOAT t){
		switch(n){
		case 2:
			return bezier_evaluate_deriv_2(p, t);
		case 4:
			return bezier_evaluate_deriv_4(p, t);
		default:
			{
				alloca_vector<FLOAT> B(n);
				bernstein_deriv_n(t,&B[0],n);
				return mul_coef(&B[0],p,n);
			}
		}
    }

#define DEC_BEZIER_LINE_EVALUATE_DERV(TYPE,FLOAT)                                               \
  TYPE bezier_evaluate_deriv(const TYPE p[4], FLOAT t){return bezier_evaluate_deriv_4(p,t);}\
  TYPE bezier_evaluate_deriv(const TYPE p[4], const FLOAT coe[4]){return bezier_evaluate_deriv_4(p,coe);}\
  TYPE bezier_evaluate_deriv(const TYPE p[] ,int n, FLOAT t){return bezier_evaluate_deriv_n(p,n,t);}

  DEC_BEZIER_LINE_EVALUATE_DERV(float,float)
  //DEC_BEZIER_LINE_EVALUATE_DERV(vector2f,float)
  DEC_BEZIER_LINE_EVALUATE_DERV(vector3f,float)
  //DEC_BEZIER_LINE_EVALUATE_DERV(vector4f,float)
  DEC_BEZIER_LINE_EVALUATE_DERV(double,double)
  //DEC_BEZIER_LINE_EVALUATE_DERV(vector2d,double)
  DEC_BEZIER_LINE_EVALUATE_DERV(vector3d,double)
  //DEC_BEZIER_LINE_EVALUATE_DERV(vector4d,double)

#undef DEC_BEZIER_LINE_EVALUATE_DERV

    //-----------------------------------------------
    template<class T, class FLOAT>
    static inline T bezier_rate(const T& a, const T& b, FLOAT t)
    {
        return (1-t)*a+t*b;
    }
	
	template<class T, class FLOAT>
    static inline void bezier_split_2(T a[2], T b[2], const T p[2], FLOAT t)
    {
        T p10 = bezier_rate(p[0],p[1],t);

        a[0] = p[0];
        a[1] = p10;

        b[0] = p10;
        b[1] = p[1];
    }
	
	template<class T, class FLOAT>
    static inline void bezier_split_3(T a[3], T b[3], const T p[3], FLOAT t)
    {
        T p10 = bezier_rate(p[0],p[1],t);
        T p11 = bezier_rate(p[1],p[2],t);

        T p20 = bezier_rate(p10,p11,t);

        a[0] = p[0];
        a[1] = p10;
        a[2] = p20;

        b[0] = p20;
        b[1] = p11;
        b[2] = p[2];
    }
  
    template<class T, class FLOAT>
    static inline void bezier_split_4(T a[4], T b[4], const T p[4], FLOAT t)
    {
        T p10 = bezier_rate(p[0],p[1],t);
        T p11 = bezier_rate(p[1],p[2],t);
        T p12 = bezier_rate(p[2],p[3],t);

        T p20 = bezier_rate(p10,p11,t);
        T p21 = bezier_rate(p11,p12,t);

        T p30 = bezier_rate(p20,p21,t);

        a[0] = p[0];
        a[1] = p10;
        a[2] = p20;
        a[3] = p30;

        b[0] = p30;
        b[1] = p21;
        b[2] = p12;
        b[3] = p[3];
    }

    template<class T, class FLOAT>
    static inline void bezier_split_n(T a[], T b[], const T p[], int i, int n, FLOAT t)
    {
        int tn = n-1-i;//3

        a[i] = p[0];
        b[tn] = p[tn];
        switch(tn){
        case 0:
          return;
        case 1:
        case 2:
        case 3:
       	case 4:
          {
            T tmp[4];
            //---------------------
            for(int j = 0;j<tn;j++){
              tmp[j] = (1-t)*p[j] + t*p[j+1];//p[j] + t*(p[j+1] - p[j])
            }
            //---------------------
            bezier_split_n(a, b, &tmp[0], i+1, n, t);
          }
          break;
        default:
          {
            alloca_vector<T> tmp(tn);
            //---------------------
            for(int j = 0;j<tn;j++){
                tmp[j] = (1-t)*p[j] + t*p[j+1];//p[j] + t*(p[j+1] - p[j])
            }
            //---------------------
            bezier_split_n(a, b, &tmp[0], i+1, n, t);
          }
          break;
        }    
    }

    template<class T, class FLOAT>
    static inline void bezier_split_n(T a[], T b[],  const T p[], int n, FLOAT t)
    {
    	switch(n){
    	case 2:
    		bezier_split_2(a, b, p, t);break;
    	case 3:
    		bezier_split_3(a, b, p, t);break;
    	case 4:
    		bezier_split_4(a, b, p, t);break;
    	default:
        	bezier_split_n(a,b, p, 0, n, t);break;
    	}
    }

#ifdef _MSC_VER
#pragma warning(push)
#pragma warning(disable:4996)
#endif
	template<class T>
	static inline void transpose_matrix_4(T m[])
	{
		for(int y = 0;y<4;y++){
			for(int x = y+1;x<4;x++){
				std::swap(m[x*4+y],m[y*4+x]);
			}
		}
	}

	template<class T>
	static inline void transpose_matrix_n(T m[], int w, int h)
	{
		if(w==h){
			for(int y = 0;y<h;y++){
				for(int x = y+1;x<w;x++){
					std::swap(m[x*h+y],m[y*w+x]);
				}
			}
		}else{
			alloca_vector<T> tmp(w*h);
			std::copy(m,m+w*h,&tmp[0]);
			for(int y = 0;y<h;y++){
				for(int x = 0;x<w;x++){
					m[x*h+y] = tmp[y*w+x];
				}
			}
		}
	}
#ifdef _MSC_VER
#pragma warning(pop)
#endif

	template<class T, class FLOAT>
    static inline void bezier_split_4(T a[4*4], T b[4*4], T c[4*4], T d[4*4], const T p[4*4], FLOAT u, FLOAT v)
    {
		T tmp0[4*4];
		T tmp1[4*4];
		//split U
		for(int i = 0;i<4;i++){
			bezier_split_4(&tmp0[i*4+0], &tmp1[i*4+0], &p[i*4+0], u);
		}

		//nuuv->nvnu
		transpose_matrix_4(&tmp0[0]);
		transpose_matrix_4(&tmp1[0]);
		//

		//split V
		for(int i = 0;i<4;i++){
			bezier_split_4(&a[i*4+0], &c[i*4+0], &tmp0[i*4+0], v);
			bezier_split_4(&b[i*4+0], &d[i*4+0], &tmp1[i*4+0], v);
		}

		//nvnu->nunv
		transpose_matrix_4(a);//00
		transpose_matrix_4(b);//10
		transpose_matrix_4(c);//01
		transpose_matrix_4(d);//11
		//		
    }

	template<class T,class FLOAT>
    static inline void bezier_split_n(T a[], T b[], T c[], T d[], const T p[], int nu, int nv, FLOAT u, FLOAT v)
    {
    	if(nu==4&&nv==4){
			bezier_split_4(a, b, c, d, p, u, v);
		}else{
			alloca_vector<T> tmp0(nu*nv);//TODO
			alloca_vector<T> tmp1(nu*nv);//TODO
			//split U
			for(int i = 0;i<nv;i++){
				bezier_split_n(&tmp0[i*nu+0], &tmp1[i*nu+0], &p[i*nu+0], nu, u);
			}
	
			//nuuv->nvnu
			transpose_matrix_n(&tmp0[0],nu,nv);// 00 01
			transpose_matrix_n(&tmp1[0],nu,nv);// 10 11
			//
	
			//split V
			for(int i = 0;i<nu;i++){
				bezier_split_n(&a[i*nv+0], &c[i*nv+0], &tmp0[i*nv+0], nv, v);
				bezier_split_n(&b[i*nv+0], &d[i*nv+0], &tmp1[i*nv+0], nv, v);
			}
	
			//nvnu->nunv
			transpose_matrix_n(a,nv,nu);//00
			transpose_matrix_n(b,nv,nu);//10
			transpose_matrix_n(c,nv,nu);//01
			transpose_matrix_n(d,nv,nu);//11
			//
		}
    }

	template<class T,class FLOAT>
	static inline void bezier_split_u_n(T a[], T b[], const T p[], int nu, int nv, FLOAT u)
	{
		for(int i = 0;i<nv;i++){
			bezier_split_n(&a[i*nu+0], &b[i*nu+0], &p[i*nu+0], nu, u);
		}
	}

	template<class T,class FLOAT>
	static inline void bezier_split_v_n(T a[], T b[], const T p[], int nu, int nv, FLOAT v)
	{
        if(nv<=4){
            T tmp[4];
            T atmp[4];
            T btmp[4];
            for(int i = 0;i<nu;i++){
                for(int j = 0;j<nv;j++){
                    tmp[j] =  p[j*nu+i];
                }
                bezier_split_n(&atmp[0], &btmp[0], &tmp[0], nv, v);
                for(int j = 0;j<nv;j++){
                    a[j*nu+i] = atmp[j];
                    b[j*nu+i] = btmp[j];
                }
            }
        }else if(nv<=16){
            T tmp[16];
            T atmp[16];
            T btmp[16];
            for(int i = 0;i<nu;i++){
                for(int j = 0;j<nv;j++){
                    tmp[j] =  p[j*nu+i];
                }
                bezier_split_n(&atmp[0], &btmp[0], &tmp[0], nv, v);
                for(int j = 0;j<nv;j++){
                    a[j*nu+i] = atmp[j];
                    b[j*nu+i] = btmp[j];
                }
            }
        }else{
            alloca_vector<T> tmp(nv);
            alloca_vector<T> atmp(nv);
            alloca_vector<T> btmp(nv);
            for(int i = 0;i<nu;i++){
                for(int j = 0;j<nv;j++){
                    tmp[j] =  p[j*nu+i];
                }
                bezier_split_n(&atmp[0], &btmp[0], &tmp[0], nv, v);
                for(int j = 0;j<nv;j++){
                    a[j*nu+i] = atmp[j];
                    b[j*nu+i] = btmp[j];
                }
            }
        }
	}

    template<class T,class FLOAT>
    static inline void bezier_crop_4(T a[], const T p[], FLOAT t0, FLOAT t1)
    {
		T tmp0[4];
		T tmp1[4];
		bezier_split_4(&tmp0[0], &tmp1[0], p, t1);
		bezier_split_4(&tmp1[0], a, &tmp0[0], t0/t1);
	}

	template<class T,class FLOAT>
    static inline void bezier_crop_n(T a[], const T p[], int n, FLOAT t0, FLOAT t1)
    {
    	if(n==4){
    		bezier_crop_4(a, p, t0, t1);
    	}else if(n<=16){
            T tmp0[16];
		    T tmp1[16];
		    bezier_split_n(&tmp0[0], &tmp1[0], p, n, t1);
		    bezier_split_n(&tmp1[0], a, &tmp0[0],  n, t0/t1);
        }else{
		    alloca_vector<T> tmp0(n);
		    alloca_vector<T> tmp1(n);
		    bezier_split_n(&tmp0[0], &tmp1[0], p, n, t1);
		    bezier_split_n(&tmp1[0], a, &tmp0[0],  n, t0/t1);
    	}
	}

	template<class T, class FLOAT>
	static inline void bezier_crop_u_n(T a[], const T p[], int nu, int nv, FLOAT u0, FLOAT u1)
	{
        if(nu==4){
            for(int i = 0;i<nv;i++){
			    bezier_crop_4(&a[i*nu+0], &p[i*nu+0], u0, u1);
		    }
        }else{
		    for(int i = 0;i<nv;i++){
			    bezier_crop_n(&a[i*nu+0], &p[i*nu+0], nu, u0, u1);
		    }
        }
	}

	template<class T, class FLOAT>
	static inline void bezier_crop_v_n(T a[], const T p[], int nu, int nv, FLOAT v0, FLOAT v1)
	{
        if(nv==4){
            T tmp[4];
			T out[4];
            for(int i = 0;i<nu;i++){    
			    for(int j = 0;j<nv;j++)tmp[j]=p[j*nu+i];
			    bezier_crop_4(&out[0], &tmp[0], v0, v1);
			    for(int j = 0;j<nv;j++)a[j*nu+i]=out[j];
		    }
        }else if(nv<=16){
            T tmp[16];
		    T out[16];
            for(int i = 0;i<nu;i++){    
			    for(int j = 0;j<nv;j++)tmp[j]=p[j*nu+i];
			    bezier_crop_n(&out[0], &tmp[0], nv, v0, v1);
			    for(int j = 0;j<nv;j++)a[j*nu+i]=out[j];
		    }
        }else{
            alloca_vector<T> tmp(nv);
			alloca_vector<T> out(nv);
		    for(int i = 0;i<nu;i++){    
			    for(int j = 0;j<nv;j++)tmp[j]=p[j*nu+i];
			    bezier_crop_n(&out[0], &tmp[0], nv, v0, v1);
			    for(int j = 0;j<nv;j++)a[j*nu+i]=out[j];
		    }
        }
	}



//--------------------------------------------------------------------------------

#define DEC_BEZIER_SPLIT(TYPE,FLOAT)                                                    \
    void bezier_split(TYPE a[4], TYPE b[4], const TYPE p[4], FLOAT t){             \
        bezier_split_4(a, b, p, t);                                               \
    }                                                                             \
    void bezier_split(TYPE a[], TYPE b[],  const TYPE p[], int n, FLOAT t){        \
        bezier_split_n(a,b,p,n,t);                                                \
    }                                                                             \
	void bezier_split(TYPE a[4*4], TYPE b[4*4], TYPE c[4*4], TYPE d[4*4], const TYPE p[4*4],  FLOAT u, FLOAT v){      \
		bezier_split_4(a,b,c,d,p,u,v);                                                                              \
	}                                                                                                               \
	void bezier_split(TYPE a[4][4], TYPE b[4][4], TYPE c[4][4], TYPE d[4][4], const TYPE p[4][4],  FLOAT u, FLOAT v){ \
		bezier_split_4((TYPE*)a,(TYPE*)b,(TYPE*)c,(TYPE*)d,(const TYPE*)p, u, v);                                   \
	}                                                                                                          \
	void bezier_split(TYPE a[], TYPE b[], TYPE c[], TYPE d[], const TYPE p[], int nu, int nv, FLOAT u, FLOAT v){ \
        bezier_split_n(a,b,c,d,p,nu,nv,u,v);                                                                   \
    }                                                                                                          \
	void bezier_split_u(TYPE a[], TYPE b[], const TYPE p[], int nu, int nv, FLOAT u){ \
		bezier_split_u_n(a, b, p, nu, nv, u);                                        \
	}                                                                                \
	void bezier_split_v(TYPE a[], TYPE b[], const TYPE p[], int nu, int nv, FLOAT v){ \
		bezier_split_v_n(a, b, p, nu, nv, v);                                        \
	}
  
  DEC_BEZIER_SPLIT(float,float)
  //DEC_BEZIER_SPLIT(vector2f,float)
  DEC_BEZIER_SPLIT(vector3f,float)
  //DEC_BEZIER_SPLIT(vector4f,float)
  DEC_BEZIER_SPLIT(double,double)
  //DEC_BEZIER_SPLIT(vector2d,double)
  DEC_BEZIER_SPLIT(vector3d,double)
  //DEC_BEZIER_SPLIT(vector4d,double)

#undef DEC_BEZIER_SPLIT

#define DEC_BEZIER_CROP(TYPE,FLOAT)                                                        \
    void bezier_crop_u(TYPE a[], const TYPE p[], int nu, int nv, FLOAT u0, FLOAT u1){  \
        bezier_crop_u_n(a, p, nu, nv, u0, u1);									     \
	}                                                                                \
	void bezier_crop_v(TYPE a[], const TYPE p[], int nu, int nv, FLOAT v0, FLOAT v1){  \
        bezier_crop_v_n(a, p, nu, nv, v0, v1);									     \
	}

  DEC_BEZIER_CROP(float,float)
  //DEC_BEZIER_CROP(vector2f,float)
  DEC_BEZIER_CROP(vector3f,float)
  //DEC_BEZIER_CROP(vector4f,float)
  DEC_BEZIER_CROP(double,double)
  //DEC_BEZIER_CROP(vector2d,double)
  DEC_BEZIER_CROP(vector3d,double)
  //DEC_BEZIER_CROP(vector4d,double)

#undef DEC_BEZIER_CROP


//--------------------------------------------------------------------------------
  template<class T, class FLOAT>
  static inline T bezier_evaluate_4(const T p[4*4],  FLOAT u, FLOAT v){
        FLOAT Bu[4];
        FLOAT Bv[4];
        bernstein(u,Bu);
        bernstein(v,Bv);

        return 
            mul_coef(
                Bv,
                mul_coef(Bu,p[ 0],p[ 1],p[ 2],p[ 3]),
                mul_coef(Bu,p[ 4],p[ 5],p[ 6],p[ 7]),
                mul_coef(Bu,p[ 8],p[ 9],p[10],p[11]),
                mul_coef(Bu,p[12],p[13],p[14],p[15])
            );    
  }
  template<class T, class FLOAT>
  static inline T bezier_evaluate_n(const T p[], int nu, int nv, FLOAT u, FLOAT v){
      if(nu==4&&nv==4){
          return bezier_evaluate_4(p, u, v);
      }else if(nu<=16&&nv<=16){
          FLOAT Bu[16];
          FLOAT Bv[16];
          bernstein(u,&Bu[0],nu);
          bernstein(v,&Bv[0],nv);

          T pv[16];
          for(int j=0;j<nv;j++){
              pv[j] = mul_coef(&Bu[0], &p[j*nu], nu);
          }

          return mul_coef(&Bv[0],&pv[0],nv);
      }else{
          alloca_vector<FLOAT> Bu(nu);
          alloca_vector<FLOAT> Bv(nv);
          bernstein(u,&Bu[0],nu);
          bernstein(v,&Bv[0],nv);

          alloca_vector<T> pv(nv);
          for(int j=0;j<nv;j++){
              pv[j] = mul_coef(&Bu[0], &p[j*nu], nu);
          }

          return mul_coef(&Bv[0],&pv[0],nv);
      }
  }

#define DEC_BEZIER_PATCH_EVALUATE(TYPE,FLOAT)                     \
  TYPE bezier_evaluate(const TYPE p[4*4],  FLOAT u, FLOAT v){ \
      return bezier_evaluate_4(p,u,v);                      \
  }                                                         \
  TYPE bezier_evaluate(const TYPE p[4][4], FLOAT u, FLOAT v){ \
      return bezier_evaluate((const TYPE*)p,u,v);           \
  }                                                         \
  TYPE bezier_evaluate(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v){\
      return bezier_evaluate_n(p,nu,nv,u,v);                \
  }
  
  DEC_BEZIER_PATCH_EVALUATE(float, float)
  //DEC_BEZIER_PATCH_EVALUATE(vector2f, float)
  DEC_BEZIER_PATCH_EVALUATE(vector3f, float)
  //DEC_BEZIER_PATCH_EVALUATE(vector4f, float)
  DEC_BEZIER_PATCH_EVALUATE(double, double)
  //DEC_BEZIER_PATCH_EVALUATE(vector2d, double)
  DEC_BEZIER_PATCH_EVALUATE(vector3d, double)
  //DEC_BEZIER_PATCH_EVALUATE(vector4d, double)

#undef DEC_BEZIER_PATCH_EVALUATE

//--------------------------------------------------------------------------------
/*
    void    bezier_evaluate_u(vector3 o[4], const vector3 p[4*4], FLOAT u[4]){
        o[0] = mul_coef(u,p[ 0],p[ 1],p[ 2],p[ 3]);
        o[1] = mul_coef(u,p[ 4],p[ 5],p[ 6],p[ 7]);
        o[2] = mul_coef(u,p[ 8],p[ 9],p[10],p[11]);
        o[3] = mul_coef(u,p[12],p[13],p[14],p[15]);
    }

    void    bezier_evaluate_v(vector3 o[4], const vector3 p[4*4], FLOAT v[4]){
        o[0] = mul_coef(v,p[ 0],p[ 4],p[ 8],p[12]);
        o[1] = mul_coef(v,p[ 1],p[ 5],p[ 9],p[13]);
        o[2] = mul_coef(v,p[ 2],p[ 6],p[10],p[14]);
        o[3] = mul_coef(v,p[ 3],p[ 7],p[11],p[15]);
    }
*/
//--------------------------------------------------------------------------------
    template<class T, class FLOAT>
    static inline T bezier_evaluate_deriv_u_4(const T p[4*4], FLOAT u, FLOAT v){
        FLOAT Bu[4];
        FLOAT Bv[4];
        bernstein_deriv_4(u,Bu);
        bernstein_4      (v,Bv);

        return mul_coef(
                Bu,
                mul_coef(Bv,p[ 0],p[ 4],p[ 8],p[12]),
                mul_coef(Bv,p[ 1],p[ 5],p[ 9],p[13]),
                mul_coef(Bv,p[ 2],p[ 6],p[10],p[14]),
                mul_coef(Bv,p[ 3],p[ 7],p[11],p[15])
            );
    }
    
  template<class T, class FLOAT>
  static inline T bezier_evaluate_deriv_v_4(const T p[4*4], FLOAT u, FLOAT v){
        FLOAT Bu[4];
        FLOAT Bv[4];
        bernstein_4(u,Bu);
        bernstein_deriv_4(v,Bv);

        return mul_coef(
            Bv,
            mul_coef(Bu,p[ 0],p[ 1],p[ 2],p[ 3]),
            mul_coef(Bu,p[ 4],p[ 5],p[ 6],p[ 7]),
            mul_coef(Bu,p[ 8],p[ 9],p[10],p[11]),
            mul_coef(Bu,p[12],p[13],p[14],p[15])
        );
  }

  template<class T, class FLOAT>
  static inline T bezier_evaluate_deriv_u_n(const T p[], int nu, int nv, FLOAT u, FLOAT v){
      if(nu==4&&nv==4){
          return bezier_evaluate_deriv_u_4(p, u, v);
      }else if(nu<=16&&nv<=16){
          FLOAT Bu[16];
          FLOAT Bv[16];
          bernstein_deriv_n(u,&Bu[0],nu);
          bernstein_n      (v,&Bv[0],nv);

          T pu[16];
          T pv[16];
          for(int j = 0;j<nu;j++){
            for(int i = 0;i<nv;i++){
              pv[i] = p[i*nu+j];
            }
            pu[j] = mul_coef(&Bv[0],&pv[0],nv);
          }

          return mul_coef(&Bu[0],&pu[0],nu);

      }else{
          alloca_vector<FLOAT> Bu(nu);
          alloca_vector<FLOAT> Bv(nv);
          bernstein_deriv_n(u,&Bu[0],nu);
          bernstein_n      (v,&Bv[0],nv);

          alloca_vector<T> pu(nu);
          alloca_vector<T> pv(nv);
          for(int j = 0;j<nu;j++){
            for(int i = 0;i<nv;i++){
              pv[i] = p[i*nu+j];
            }
            pu[j] = mul_coef(&Bv[0],&pv[0],nv);
          }

          return mul_coef(&Bu[0],&pu[0],nu);
      }
  }
    
  template<class T, class FLOAT>
  static inline T bezier_evaluate_deriv_v_n(const T p[], int nu, int nv, FLOAT u, FLOAT v){
      if(nu==4&&nv==4){
          return bezier_evaluate_deriv_v_4(p, u, v);
      }else if(nu<=16&&nv<=16){
          FLOAT Bu[16];
          FLOAT Bv[16];
          bernstein_n      (u,&Bu[0],nu);
          bernstein_deriv_n(v,&Bv[0],nv);

          T pu[16];
          T pv[16];
          for(int i = 0;i<nv;i++){
            for(int j = 0;j<nu;j++){
              pu[j] = p[i*nu+j];
            }
            pv[i] = mul_coef(&Bu[0],&pu[0],nu);
          }

          return mul_coef(&Bv[0],&pv[0],nv);
      }else{
          alloca_vector<FLOAT> Bu(nu);
          alloca_vector<FLOAT> Bv(nv);
          bernstein_n      (u,&Bu[0],nu);
          bernstein_deriv_n(v,&Bv[0],nv);

          alloca_vector<T> pu(nu);
          alloca_vector<T> pv(nv);
          for(int i = 0;i<nv;i++){
            for(int j = 0;j<nu;j++){
              pu[j] = p[i*nu+j];
            }
            pv[i] = mul_coef(&Bu[0],&pu[0],nu);
          }

          return mul_coef(&Bv[0],&pv[0],nv);
      }
  }

#define DEC_BEZIER_EVALUATE_DPDU(TYPE, FLOAT)                                \
    TYPE bezier_evaluate_deriv_u(const TYPE p[4*4], FLOAT u, FLOAT v)  {\
        return bezier_evaluate_deriv_u_4(p,u,v);                      \
    }                                                                 \
    TYPE bezier_evaluate_deriv_u(const TYPE p[4][4], FLOAT u, FLOAT v) {\
        return bezier_evaluate_deriv_u((const TYPE*)p,u,v);           \
    }                                                                 \
    TYPE bezier_evaluate_deriv_u(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v)  {\
        return bezier_evaluate_deriv_u_n(p,nu,nv,u,v);                \
    } 

    DEC_BEZIER_EVALUATE_DPDU(float, float)
    //DEC_BEZIER_EVALUATE_DPDU(vector2f, float)
    DEC_BEZIER_EVALUATE_DPDU(vector3f, float)
    //DEC_BEZIER_EVALUATE_DPDU(vector4f, float)
	DEC_BEZIER_EVALUATE_DPDU(double, double)
    //DEC_BEZIER_EVALUATE_DPDU(vector2d, double)
    DEC_BEZIER_EVALUATE_DPDU(vector3d, double)
    //DEC_BEZIER_EVALUATE_DPDU(vector4d, double)

#undef DEC_BEZIER_EVALUATE_DPDU

#define DEC_BEZIER_EVALVATE_DPDV(TYPE, FLOAT)                           \
    TYPE bezier_evaluate_deriv_v(const TYPE p[4*4], FLOAT u, FLOAT v)  {\
        return bezier_evaluate_deriv_v_4(p,u,v);                        \
    }                                                                   \
    TYPE bezier_evaluate_deriv_v(const TYPE p[4][4], FLOAT u, FLOAT v) {\
        return bezier_evaluate_deriv_v((const TYPE*)p,u,v);             \
    }                                                                   \
    TYPE bezier_evaluate_deriv_v(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v)  {\
        return bezier_evaluate_deriv_v_n(p,nu,nv,u,v);                \
    } 

    DEC_BEZIER_EVALVATE_DPDV(float, float)
    //DEC_BEZIER_EVALVATE_DPDV(vector2f, float)
    DEC_BEZIER_EVALVATE_DPDV(vector3f, float)
    //DEC_BEZIER_EVALVATE_DPDV(vector4f, float)
	DEC_BEZIER_EVALVATE_DPDV(double, double)
    //DEC_BEZIER_EVALVATE_DPDV(vector2d, double)
    DEC_BEZIER_EVALVATE_DPDV(vector3d, double)
    //DEC_BEZIER_EVALVATE_DPDV(vector4d, double)

#undef DEC_BEZIER_EVALUATE_DPDV

//--------------------------------------------------------------------------------

	template<class T, class FLOAT>
    static inline T bezier_evaluate_normal_(const T p[4*4], FLOAT u, FLOAT v)
    {
        T sU = bezier_evaluate_deriv_u_4(p,u,v);
        T sV = bezier_evaluate_deriv_v_4(p,u,v);       
        return normalize(cross(sU,sV));
    }

	template<class T, class FLOAT>
    static inline T bezier_evaluate_normal_(const T p[4][4], FLOAT u, FLOAT v){
        return bezier_evaluate_normal((const T*)p,u,v);
    }

	template<class T, class FLOAT>
    static inline T bezier_evaluate_normal_(const T p[], int nu, int nv, FLOAT u, FLOAT v)
    {
        T sU = bezier_evaluate_deriv_u_n(p,nu,nv,u,v);
        T sV = bezier_evaluate_deriv_v_n(p,nu,nv,u,v);       
        return normalize(cross(sU,sV));
    }

#define DEC_BEZIER_EVALUATE_NORMAL(TYPE,FLOAT)                         \
	TYPE bezier_evaluate_normal(const TYPE p[4*4], FLOAT u, FLOAT v){  \
		return bezier_evaluate_normal_(p, u, v);                       \
	}                                                                  \
	TYPE bezier_evaluate_normal(const TYPE p[4][4], FLOAT u, FLOAT v){ \
		return bezier_evaluate_normal_(p, u, v);                       \
	}                                                                  \
	TYPE bezier_evaluate_normal(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v){\
		return bezier_evaluate_normal_(p, nu, nv, u, v);                          \
	}

	DEC_BEZIER_EVALUATE_NORMAL(vector3f,float)
	DEC_BEZIER_EVALUATE_NORMAL(vector3d,double)

#undef DEC_BEZIER_EVALUATE_NORMAL

	template<class T>
	static inline T bezier_min_n(const T p[], int n, int d){
		T tmp = p[0];
		for(int i = 1;i<n;i++){
			for(int j = 0;j<d;j++){
				if(tmp[j]>p[i][j])tmp[j]=p[i][j];
			}
		}
		return tmp;
	}

	template<class T>
	static inline T bezier_max_n(const T p[], int n, int d){
		T tmp = p[0];
		for(int i = 1;i<n;i++){
			for(int j = 0;j<d;j++){
				if(tmp[j]<p[i][j])tmp[j]=p[i][j];
			}
		}
		return tmp;
	}

#define DEC_BEZIER_MINMAX(FLOAT)              \
	FLOAT bezier_min(const FLOAT p[], int n){ \
		FLOAT tmp = p[0];                     \
		for(int i = 1;i<n;i++){               \
			if(tmp>p[i])tmp=p[i];             \
		}                                     \
		return tmp;                           \
	}                                         \
	FLOAT bezier_max(const FLOAT p[], int n){ \
		FLOAT tmp = p[0];                     \
		for(int i = 1;i<n;i++){               \
			if(tmp<p[i])tmp=p[i];             \
		}                                     \
		return tmp;                           \
	}

	DEC_BEZIER_MINMAX(float)
	DEC_BEZIER_MINMAX(double)

#undef DEC_BEZIER_MINMAX

#define DEC_BEZIER_MINMAX(TYPE,D)             \
	TYPE bezier_min(const TYPE p[], int n){   \
		return bezier_min_n(p,n,D);           \
	}                                         \
	TYPE bezier_max(const TYPE p[], int n){   \
		return bezier_max_n(p,n,D);           \
	}

	DEC_BEZIER_MINMAX(vector2f,2)
	DEC_BEZIER_MINMAX(vector3f,3)
	//DEC_BEZIER_MINMAX(vector4f,4)
	DEC_BEZIER_MINMAX(vector2d,2)
	DEC_BEZIER_MINMAX(vector3d,3)
	//DEC_BEZIER_MINMAX(vector4d,4)

#undef DEC_BEZIER_MINMAX

}
