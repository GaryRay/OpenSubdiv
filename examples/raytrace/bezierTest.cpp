#include <cstdio>
#include <cstdlib>
#include <vector>
#include "osdutil/math.h"
#include "osdutil/math_sse.h"
#include "osdutil/bezier.h"

#include "common.h"
#include "bezier_patch.hpp"

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
        printf("(%f, %f, %f) != (%f, %f, %f), delta = (%g, %g, %g)\n",
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

template <typename PATCH, typename VEC>
static void rotateTest()
{
    typedef PATCH Patch;

    printf("Rotate Test\n");
    srand(0);

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
        //dump(patch0);
        //dump(patch1);

        Patch result0, result1;
        double min = 0.3, max = 0.6;

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

        dump(result0);
        dump(result1);

        comparePatch(result0, result1);
    }
    printf("Rotate Test end\n");
}

template <typename PATCH, typename VEC>
static void evalTest()
{
    typedef PATCH Patch;

    printf("Eval Test\n");
    srand(0);

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
    rotateTest<MallieBezierPatch, mallie::vector3>();
    evalTest<MallieBezierPatch, mallie::vector3>();

    printf("osd float\n");
    rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3f, float>, OsdUtil::vec3f>();
    evalTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3f, float>, OsdUtil::vec3f>();

    printf("osd sse\n");
    rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3sse, float>, OsdUtil::vec3sse>();
    evalTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3sse, float>, OsdUtil::vec3sse>();

    printf("osd double\n");
    rotateTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3d, double>, OsdUtil::vec3d>();
    evalTest<OsdUtil::OsdUtilBezierPatch<OsdUtil::vec3d, double>, OsdUtil::vec3d>();
}
