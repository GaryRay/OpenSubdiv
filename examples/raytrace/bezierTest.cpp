#include <cstdio>
#include <cstdlib>
#include <vector>
#include "osdutil/math.h"
#include "osdutil/math_sse.h"
#include "osdutil/bezier.h"

#include "common.h"
#include "bezier_patch.hpp"

namespace {

void
EvalCubicBSplineFloat(float u, float N[4], float NU[4])
{
  float T = u;
  float S = 1.0f - u;

  N[0] = S * S * S;
  N[1] = (4.0f * S * S * S + T * T * T) + (12.0f * S * T * S + 6.0f * T * S * T);
  N[2] = (4.0f * T * T * T + S * S * S) + (12.0f * T * S * T + 6.0f * S * T * S);
  N[3] = T * T * T;

  NU[0] = -S * S;
  NU[1] = -T * T - 4.0f * T * S;
  NU[2] =  S * S + 4.0f * S * T;
  NU[3] =  T * T;
}


// bitwise float value eq 
bool my_float_eq(float x, float y)
{
  return (*reinterpret_cast<int*>(&x)) == (*reinterpret_cast<int*>(&y));
}

#define my_eq_fassert(x, y) { \
  if (!my_float_eq((x), (y))) { \
    fprintf(stderr, "assert failed. x = %f(0x%08x), y = %f(0x%08x). line:%d\n", x, *reinterpret_cast<unsigned int *>(&(x)), y, *reinterpret_cast<unsigned int *>(&(y)), __LINE__); \
    exit(1); \
  } \
}

#define my_bool_assert(ret) { \
  if (!(ret)) { \
    fprintf(stderr, "assert failed. line: %d\n", __LINE__); \
    exit(1); \
  } \
}

void
TestCubic()
{
  // 0.0 ~ 1.0. sign = 0
  for (int e = 125; e < 126; e++) { // exponential
    printf("e = %d\n", e);
    for (int m = 0; m < 83886082; m++) { // mantissa: 8388608 = 2**23
      unsigned int ui = (e << 23) + m;
      float f = (*reinterpret_cast<float *>(&ui));
      printf("f = %f, %d\n", f, f <= 1.0f);
      my_bool_assert(f <= 1.0f);
      my_bool_assert(f >= 0.0f);

      float N0[4], N1[4];
      float NU0[4], NU1[4];
      EvalCubicBSplineFloat(f, N0, NU0);
      EvalCubicBSplineFloat(1.0f - f, N1, NU1);

      //printf("%f, %f\n", N0[0], N1[3]);
      my_eq_fassert(N0[0], N1[3]);
      my_eq_fassert(N0[1], N1[2]);
      my_eq_fassert(N0[2], N1[1]);
      my_eq_fassert(N0[3], N1[0]);

      //my_eq_fassert(NU0[0], NU1[3]);
      //my_eq_fassert(NU0[1], NU1[2]);
      //my_eq_fassert(NU0[2], NU1[1]);
      //my_eq_fassert(NU0[3], NU1[0]);
    }
  }
}


}


struct MallieBezierPatch : public mallie::bezier_patch<mallie::vector3> {
    MallieBezierPatch() :
        mallie::bezier_patch<mallie::vector3>(4, 4) {}
    MallieBezierPatch(const mallie::vector3 *p) :
        mallie::bezier_patch<mallie::vector3>(4, 4, p) { }

    mallie::vector3 Get(int i, int j) const {
        return get_cp_at(i, j);
    }
    mallie::vector3 Evaluate(float u, float v) const {
        return evaluate(u, v);
    }
    void Transpose() {
        for (int y = 0; y < get_nv(); ++y) {
            for (int x = y+1; x < get_nu(); ++x) {
                mallie::vector3 tmp = get_cp_at(x, y);
                set_cp_at(x, y, get_cp_at(y, x));
                set_cp_at(y, x, tmp);
            }
        }
    }
    void Rotate() {
        Transpose();
        for (int i = 0; i < get_nv(); ++i) {
            for (int j = 0; j < get_nu()/2; ++j) {
                mallie::vector3 tmp = get_cp_at(j, i);
                set_cp_at(j, i, get_cp_at(get_nu()-1-j, i));
                set_cp_at(get_nu()-1-j, i, tmp);
            }
        }
    }
    void CropU(MallieBezierPatch &patch, float u0, float u1) const {
        crop_u(patch, u0, u1);
    }
    void CropV(MallieBezierPatch &patch, float v0, float v1) const {
        crop_v(patch, v0, v1);
    }
};

template <typename PATCH>
static void dump(PATCH const &patch)
{
    for (int j = 0; j < 4; ++j) {
        for (int i = 0; i < 4; ++i) {
            printf("(%.8f, %.8f, %.8f) ", patch.Get(i, j)[0], patch.Get(i, j)[1], patch.Get(i, j)[2]);
        }
        printf("\n");
    }
    printf("\n");
}

template <typename VEC>
static void compare(VEC const &v0, VEC const &v1)
{
    bool fail = false;
    for (int k = 0; k < 3; ++k) {
        if (v0[k] != v1[k]) {
            fail = true;
        }
    }
    if (fail) {
        printf("(%.10f, %.10f, %.10f) != (%.10f, %.10f, %.10f), delta = (%g, %g, %g)\n",
               v0[0], v0[1], v0[2],
               v1[0], v1[1], v1[2],
               v1[0]-v0[0], v1[1]-v0[1], v1[2]-v0[2]);
    }
}

template <typename PATCH>
static void comparePatch(PATCH const &patch0, PATCH const &patch1)
{
    for (int i = 0; i < 4; ++i) {
        for (int j = 0; j < 4; ++j) {
            compare(patch0.Get(i, j), patch1.Get(i, j));
        }
    }
}

template <typename PATCH, typename VEC, typename REAL>
static void rotateTest(int seed=0)
{
    typedef PATCH Patch;

    printf("Rotate Test\n");
    srand(seed);

    VEC cp[16];
    for (int i = 0; i < 16; ++i) {
        double x = std::rand()/double(RAND_MAX);
        double y = std::rand()/double(RAND_MAX);
        double z = std::rand()/double(RAND_MAX);
        cp[i] = VEC(x, y, z);
    }

    for (int rotate = 1; rotate < 4; ++rotate) {
        Patch patch0(cp);
        Patch patch1(cp);

        for (int i = 0; i < rotate; ++i) {
            // rotate 90, 180, 270
            patch1.Rotate();
        }
        printf("Rotate = %d\n", rotate);
        // dump(patch0);
        // dump(patch1);

        Patch result0, result1;
        REAL min = 0.3, max = 0.6;

        patch0.CropU(result0, min, max);
        if (rotate == 1) {
            patch1.CropV(result1, min, max);
        } else if (rotate == 2) {
            patch1.CropU(result1, 1-max, 1-min);
        } else {
            patch1.CropV(result1, 1-max, 1-min);
        }

        for (int i = 0; i < 4-rotate; ++i) {
            // rotate 270, 180, 90
            result1.Rotate();
        }
        // dump(result0);
        // dump(result1);

        comparePatch(result0, result1);
    }
    printf("Rotate Test end\n");
}

template <typename PATCH, typename VEC>
static void evalTest(int seed = 0)
{
    typedef PATCH Patch;

    printf("Eval Test\n");
    srand(seed);

    VEC cp[16];
    for (int i = 0; i < 16; ++i) {
        double x = std::rand()/double(RAND_MAX);
        double y = std::rand()/double(RAND_MAX);
        double z = std::rand()/double(RAND_MAX);
        cp[i] = VEC(x, y, z);
    }

    for (int rotate = 1; rotate < 4; ++rotate) {
        Patch patch0(cp);
        Patch patch1(cp);

        for (int i = 0; i < rotate; ++i) {
            // rotate 90, 180, 270
            patch1.Rotate();
        }

        for (int j = 0; j < 10; ++j) {
            double u0 = std::rand()/double(RAND_MAX);
            double v0 = std::rand()/double(RAND_MAX);
            double u = u0;
            double v = v0;

            if (rotate == 1) {
                std::swap(u, v);
                u = 1.0 - u;
            } else if (rotate == 2) {
                u = 1.0 - u;
                v = 1.0 - v;
            } else if (rotate == 3) {
                std::swap(u, v);
                v = 1.0 - v;
            }

            VEC p0 = patch0.Evaluate(u0, v0);
            VEC p1 = patch1.Evaluate(u, v);

            compare(p0, p1);
        }
    }
    printf("Eval Test end\n");
}

int main()
{
    printf("original\n");
    rotateTest<MallieBezierPatch, mallie::vector3, float>();
    //    evalTest<MallieBezierPatch, mallie::vector3>();

    printf("osd float\n");
    rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3f, float>, OsdUtil::vec3f, float>();
    //    evalTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3f, float>, OsdUtil::vec3f>();

    printf("osd sse\n");
    rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3f, float>, OsdUtil::vec3f, float>();
    //rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3sse, float>, OsdUtil::vec3sse>();
    //    evalTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3sse, float>, OsdUtil::vec3sse>();

    printf("osd double\n");
    rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3d, double>, OsdUtil::vec3d, double>();
    //    evalTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3d, double>, OsdUtil::vec3d>();

    // May take time
    // todo
    //TestCubic();
}
