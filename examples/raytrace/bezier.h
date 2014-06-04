#ifndef __BEZIER_H__
#define __BEZIER_H__

#include "vector2.h"
#include "vector3.h"

#include <stdlib.h>
#include <string.h>

namespace mallie{

#define DEF_BEZIER_BERNSTEIN(FLOAT)                 \
	void bernstein(FLOAT t, FLOAT e[4]);            \
	void bernstein(FLOAT t, FLOAT e[], int n);      \
	void bernstein_deriv(FLOAT t, FLOAT e[4]);      \
	void bernstein_deriv(FLOAT t, FLOAT e[], int n);

	DEF_BEZIER_BERNSTEIN(float)
	DEF_BEZIER_BERNSTEIN(double)

#undef DEC_BEZIER_BERNSTEIN

#define DEF_BEZIER_EVALUATE_NORMAL(TYPE,FLOAT)                         \
	TYPE bezier_evaluate_normal(const TYPE p[4*4], FLOAT u, FLOAT v);  \
	TYPE bezier_evaluate_normal(const TYPE p[4][4], FLOAT u, FLOAT v); \
	TYPE bezier_evaluate_normal(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v);

	DEF_BEZIER_EVALUATE_NORMAL(vector3f,float)
	DEF_BEZIER_EVALUATE_NORMAL(vector3d,double)

#undef DEF_BEZIER_EVALUATE_NORMAL


#define DEF_BEZIER_LINE_EVALUATE(TYPE, FLOAT)               \
  TYPE bezier_evaluate(const TYPE p[4], FLOAT t);           \
  TYPE bezier_evaluate(const TYPE p[4], const FLOAT coe[4]);\
  TYPE bezier_evaluate(const TYPE p[] ,int n, FLOAT t);

  DEF_BEZIER_LINE_EVALUATE(float,float)
  //DEF_BEZIER_LINE_EVALUATE(vector2f,float)
  DEF_BEZIER_LINE_EVALUATE(vector3f,float)
  //DEF_BEZIER_LINE_EVALUATE(vector4f,float)
  DEF_BEZIER_LINE_EVALUATE(double,double)
  //DEF_BEZIER_LINE_EVALUATE(vector2d,double)
  DEF_BEZIER_LINE_EVALUATE(vector3d,double)
  //DEF_BEZIER_LINE_EVALUATE(vector4d,double)

#undef DEF_BEZIER_LINE_EVALUATE

#define DEF_BEZIER_LINE_EVALUATE_DERV(TYPE,FLOAT)                                                              \
  TYPE bezier_evaluate_deriv(const TYPE p[4], FLOAT t);                                                        \
  TYPE bezier_evaluate_deriv(const TYPE p[4], const FLOAT coe[4]);                                             \
  TYPE bezier_evaluate_deriv(const TYPE p[] ,int n, FLOAT t);                                                  \
  inline TYPE bezier_evaluate_dPdT(const TYPE p[4], FLOAT t){return bezier_evaluate_deriv(p, t);}              \
  inline TYPE bezier_evaluate_dPdT(const TYPE p[4], const FLOAT coe[4]){return bezier_evaluate_deriv(p, coe);} \
  inline TYPE bezier_evaluate_dPdT(const TYPE p[] ,int n, FLOAT t){return bezier_evaluate_deriv(p, n, t);} 

  DEF_BEZIER_LINE_EVALUATE_DERV(float,float)
  //DEF_BEZIER_LINE_EVALUATE_DERV(vector2f,float)
  DEF_BEZIER_LINE_EVALUATE_DERV(vector3f,float)
  //DEF_BEZIER_LINE_EVALUATE_DERV(vector4f,float)
  DEF_BEZIER_LINE_EVALUATE_DERV(double,double)
  //DEF_BEZIER_LINE_EVALUATE_DERV(vector2d,double)
  DEF_BEZIER_LINE_EVALUATE_DERV(vector3d,double)
  //DEF_BEZIER_LINE_EVALUATE_DERV(vector4d,double)

#undef DEF_BEZIER_LINE_EVALUATE_DERV

#define DEF_BEZIER_PATCH_EVALUATE(TYPE, FLOAT)               \
  TYPE bezier_evaluate(const TYPE p[4*4],  FLOAT u, FLOAT v);\
  TYPE bezier_evaluate(const TYPE p[4][4], FLOAT u, FLOAT v);\
  TYPE bezier_evaluate(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v);

  DEF_BEZIER_PATCH_EVALUATE(float,float)
  //DEF_BEZIER_PATCH_EVALUATE(vector2f,float)
  DEF_BEZIER_PATCH_EVALUATE(vector3f,float)
  //DEF_BEZIER_PATCH_EVALUATE(vector4f,float)

  DEF_BEZIER_PATCH_EVALUATE(double,double)
  //DEF_BEZIER_PATCH_EVALUATE(vector2d,double)
  DEF_BEZIER_PATCH_EVALUATE(vector3d,double)
  //DEF_BEZIER_PATCH_EVALUATE(vector4d,double)

#undef DEF_BEZIER_PATCH_EVALUATE

#define DEF_BEZIER_PATCH_EVALUATE_DERV(TYPE,FLOAT)                                 \
	TYPE bezier_evaluate_deriv_u(const TYPE p[4*4],  FLOAT u, FLOAT v);            \
	TYPE bezier_evaluate_deriv_u(const TYPE p[4][4], FLOAT u, FLOAT v);            \
	TYPE bezier_evaluate_deriv_u(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v);\
	TYPE bezier_evaluate_deriv_v(const TYPE p[4*4],  FLOAT u, FLOAT v);            \
	TYPE bezier_evaluate_deriv_v(const TYPE p[4][4], FLOAT u, FLOAT v);            \
	TYPE bezier_evaluate_deriv_v(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v);\
    inline TYPE bezier_evaluate_dPdU(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v){return bezier_evaluate_deriv_u(p, nu, nv, u, v);}\
    inline TYPE bezier_evaluate_dPdV(const TYPE p[], int nu, int nv, FLOAT u, FLOAT v){return bezier_evaluate_deriv_v(p, nu, nv, u, v);}

	DEF_BEZIER_PATCH_EVALUATE_DERV(float, float)
	//DEF_BEZIER_PATCH_EVALUATE_DERV(vector2f,float)
	DEF_BEZIER_PATCH_EVALUATE_DERV(vector3f,float)
	//DEF_BEZIER_PATCH_EVALUATE_DERV(vector4f,float)
	DEF_BEZIER_PATCH_EVALUATE_DERV(double, double)
	//DEF_BEZIER_PATCH_EVALUATE_DERV(vector2d, double)
	DEF_BEZIER_PATCH_EVALUATE_DERV(vector3d, double)
	//DEF_BEZIER_PATCH_EVALUATE_DERV(vector4d, double)

#undef DEF_BEZIER_PATCH_EVALUATE_DERV

	//de Casteljau
#define DEF_BEZIER_SPLIT(TYPE,FLOAT)                                                        \
  void bezier_split(TYPE a[4], TYPE b[4], const TYPE p[4],        FLOAT t = FLOAT(0.5));\
  void bezier_split(TYPE a[],  TYPE b[],  const TYPE p[],  int n, FLOAT t = FLOAT(0.5));\
  void bezier_split(TYPE a[4*4], TYPE b[4*4], TYPE c[4*4], TYPE d[4*4], const TYPE p[4*4],  FLOAT u = FLOAT(0.5), FLOAT v = FLOAT(0.5));     \
  void bezier_split(TYPE a[4][4], TYPE b[4][4], TYPE c[4][4], TYPE d[4][4], const TYPE p[4][4], FLOAT u = FLOAT(0.5), FLOAT v = FLOAT(0.5)); \
  void bezier_split(TYPE a[], TYPE b[], TYPE c[], TYPE d[], const TYPE p[], int nu, int nv, FLOAT u = FLOAT(0.5), FLOAT v = FLOAT(0.5));     \
  void bezier_split_u(TYPE a[], TYPE b[], const TYPE p[], int nu, int nv, FLOAT u = FLOAT(0.5)); \
  void bezier_split_v(TYPE a[], TYPE b[], const TYPE p[], int nu, int nv, FLOAT v = FLOAT(0.5));

  DEF_BEZIER_SPLIT(float,float)
  //DEF_BEZIER_SPLIT(vector2f,float)
  DEF_BEZIER_SPLIT(vector3f,float)
  //DEF_BEZIER_SPLIT(vector4f,float)
  DEF_BEZIER_SPLIT(double, double)
  //DEF_BEZIER_SPLIT(vector2d, double)
  DEF_BEZIER_SPLIT(vector3d, double)
  //DEF_BEZIER_SPLIT(vector4d, double)


#undef DEF_BEZIER_SPLIT

#define DEF_BEZIER_CROP(TYPE,FLOAT)                                                        \
    void bezier_crop_u(TYPE a[], const TYPE p[], int nu, int nv, FLOAT u0, FLOAT u1);  \
	void bezier_crop_v(TYPE a[], const TYPE p[], int nu, int nv, FLOAT u0, FLOAT u1);

	DEF_BEZIER_CROP(float,float)
	//DEF_BEZIER_CROP(vector2f,float)
	DEF_BEZIER_CROP(vector3f,float)
	//DEF_BEZIER_CROP(vector4f,float)
	DEF_BEZIER_CROP(double,double)
	//DEF_BEZIER_CROP(vector2d,double)
	DEF_BEZIER_CROP(vector3d,double)
	//DEF_BEZIER_CROP(vector4d,double)

#undef DEC_BEZIER_CROP
 


#define DEF_BEZIER_MINMAX(TYPE)             \
	TYPE bezier_min(const TYPE a[], int n); \
	TYPE bezier_max(const TYPE a[], int n); \
  void bezier_minmax(TYPE& min, TYPE& max, const TYPE p[], int n);

	DEF_BEZIER_MINMAX(float)
	DEF_BEZIER_MINMAX(vector2f)
	DEF_BEZIER_MINMAX(vector3f)
	//DEF_BEZIER_MINMAX(vector4f)
	DEF_BEZIER_MINMAX(double)
	DEF_BEZIER_MINMAX(vector2d)
	DEF_BEZIER_MINMAX(vector3d)
	//DEF_BEZIER_MINMAX(vector4d)

#undef DEF_BEZIER_MINMAX
    
}


#endif