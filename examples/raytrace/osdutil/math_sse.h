#ifndef OSDUTIL_MATH_SSE_H
#define OSDUTIL_MATH_SSE_H

//#include "../../../opensubdiv/version.h"  // FIXME
//namespace OpenSubdiv {
//namespace OPENSUBDIV_VERSION {
#include <xmmintrin.h>

namespace OsdUtil {

struct matrix3sse;
struct matrix4sse;

#ifdef _MSC_VER
    typedef _declspec(align(16)) float float_array[4];
#else
    typedef  __attribute__ ((aligned (16))) float float_array[4];
#endif

struct vec3sse {
    typedef float ElementType;
    typedef matrix3sse Matrix3Type;
    typedef matrix4sse Matrix4Type;

    vec3sse() {
    }
    vec3sse(__m128 vec4) {
        v = vec4;
    }
    vec3sse(float f) {
        float_array mm = { f, f, f, 1};
        v = _mm_load_ps(mm);
    }
    vec3sse(float x, float y, float z) {
        float_array mm = { x, y, z, 1};
        v = _mm_load_ps(mm);
    }

    vec3sse operator*(float f) const {
        __m128 vf = { f, f, f, 1 };
        return vec3sse(_mm_mul_ps(v, vf));
    }
    vec3sse operator+(vec3sse f) const {
        return vec3sse(_mm_add_ps(v, f.v));
    }
    vec3sse operator-(vec3sse f) const {
        return vec3sse(_mm_sub_ps(v, f.v));
    }
    vec3sse &operator+=(const vec3sse &f) {
        v = _mm_add_ps(v, f.v);
        return (*this);
    }
    vec3sse &operator-=(const vec3sse &f) {
        v = _mm_sub_ps(v, f.v);
        return (*this);
    }
    float operator[](int i) const {
        float_array mm;
        _mm_store_ps(mm, v);
        return mm[i];
    }

    template <typename T>
    vec3sse &operator=(const T &s) {
        v = _mm_loadu_ps((const float*)s.v);
        // set w=1. sad.
        v = _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 1, 2, 3));
        v = _mm_move_ss(v, _mm_set_ss(1));
        v = _mm_shuffle_ps(v, v, _MM_SHUFFLE(0, 1, 2, 3));
        return *this;
    }

    float length() const {
        float_array mm;
        _mm_store_ps(mm, _mm_mul_ps(v, v));
        return std::sqrt(mm[0] + mm[1] + mm[2]);
    }

    void normalize() {
        float il = 1.0 / length();
        float_array mm = { il, il, il, 1 };
        v = _mm_mul_ps(v, _mm_load_ps(mm));
    }

    vec3sse min(vec3sse const &a) const {
        return vec3sse(_mm_min_ps(v, a.v));
    }
    vec3sse max(vec3sse const &a) const {
        return vec3sse(_mm_max_ps(v, a.v));
    }

    static vec3sse min(vec3sse const &a, vec3sse const &b) {
        return vec3sse(_mm_min_ps(a.v, b.v));
    }

    static vec3sse max(vec3sse const &a, vec3sse const &b) {
        return vec3sse(_mm_max_ps(a.v, b.v));
    }

    __m128 v;
};

inline static vec3sse cross(vec3sse const &a, vec3sse const &b) {

    __m128 r = _mm_sub_ps(
        _mm_mul_ps(_mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(3, 0, 2, 1)),
                   _mm_shuffle_ps(b.v, b.v, _MM_SHUFFLE(3, 1, 0, 2))),
        _mm_mul_ps(_mm_shuffle_ps(a.v, a.v, _MM_SHUFFLE(3, 1, 0, 2)),
                   _mm_shuffle_ps(b.v, b.v, _MM_SHUFFLE(3, 0, 2, 1))));
    return r;
}

inline static vec3sse operator*(float f, vec3sse const &v) {
    return v*f;
}

// ----------------------------------------------------------------------------

struct matrix3sse {
public:
    matrix3sse() { }
    matrix3sse(vec3sse m0, vec3sse m1, vec3sse m2) {
        v[0] = m0.v;
        v[1] = m1.v;
        v[2] = m2.v;
        for (int i = 0; i < 3; ++i) {
            v[i] = _mm_shuffle_ps(v[i], v[i], _MM_SHUFFLE(0, 1, 2, 3));
            v[i] = _mm_move_ss(v[i], _mm_setzero_ps());
            v[i] = _mm_shuffle_ps(v[i], v[i], _MM_SHUFFLE(0, 1, 2, 3));
        }
    }
    matrix3sse(float m00, float m01, float m02,
               float m10, float m11, float m12,
               float m20, float m21, float m22) {

        float_array m0 = { m00, m01, m02, 0 };
        float_array m1 = { m10, m11, m12, 0 };
        float_array m2 = { m20, m21, m22, 0 };
        v[0] = _mm_load_ps(m0);
        v[1] = _mm_load_ps(m1);
        v[2] = _mm_load_ps(m2);
    }

    // getRotate constructor.
    matrix3sse(vec3sse const &dx) {
        // [x, y, z, w] -> [x/sqrt(x*x+y*y), y/sqrt(x*x+y*y), 0, 0]
        __m128 t = _mm_mul_ps(dx.v, dx.v);
        __m128 n = _mm_add_ps(t, _mm_shuffle_ps(t, t, _MM_SHUFFLE(3, 2, 0, 1)));
        t = _mm_mul_ps(dx.v, _mm_rsqrt_ps(n));
        __m128 x = _mm_shuffle_ps(t, _mm_setzero_ps(), _MM_SHUFFLE(3, 2, 1, 0));

        // vec3sse y(-x[1], x[0], 0); // -1 1 0 0
        __m128 y = _mm_shuffle_ps(x, _mm_setzero_ps(), _MM_SHUFFLE(3, 2, 0, 1));
        __m128 sign = _mm_castsi128_ps(_mm_set_epi32(0, 0, 0, 0x80000000));
        y = _mm_xor_ps(y, sign);

        // z=w=0
        v[0] = _mm_shuffle_ps(x, _mm_setzero_ps(), _MM_SHUFFLE(3, 2, 1, 0));
        v[1] = _mm_shuffle_ps(y, _mm_setzero_ps(), _MM_SHUFFLE(3, 2, 1, 0));
        // v[2] = 0 0 1
        v[2] = _mm_movelh_ps(_mm_setzero_ps(), _mm_set_ss(1));
    }

    vec3sse operator*(vec3sse vec) const {
        __m128 t0 = _mm_mul_ps(v[0], vec.v);
        __m128 t1 = _mm_mul_ps(v[1], vec.v);
        __m128 t2 = _mm_mul_ps(v[2], vec.v);
        __m128 t3 = _mm_set_ss(1); // 0,0,0,1
        _MM_TRANSPOSE4_PS(t0, t1, t2, t3);
        return _mm_add_ps(_mm_add_ps(t0, t1), _mm_add_ps(t2, t3));
    }

    __m128 v[3];
};



// ----------------------------------------------------------------------------

struct matrix4sse {
public:
    matrix4sse() { }
    matrix4sse(__m128 r0, __m128 r1, __m128 r2, __m128 r3) {
        v[0] = r0; v[1] = r1; v[2] = r2; v[3] = r3;
    }
    matrix4sse(float m00, float m01, float m02, float m03,
               float m10, float m11, float m12, float m13,
               float m20, float m21, float m22, float m23,
               float m30, float m31, float m32, float m33) {
        float_array m0 = { m00, m01, m02, m03 };
        float_array m1 = { m10, m11, m12, m13 };
        float_array m2 = { m20, m21, m22, m23 };
        float_array m3 = { m30, m31, m32, m33 };
        v[0] = _mm_load_ps(m0);
        v[1] = _mm_load_ps(m1);
        v[2] = _mm_load_ps(m2);
        v[3] = _mm_load_ps(m3);
    }

    vec3sse operator*(vec3sse const & vec) const {
        __m128 t0 = _mm_mul_ps(v[0], vec.v);
        __m128 t1 = _mm_mul_ps(v[1], vec.v);
        __m128 t2 = _mm_mul_ps(v[2], vec.v);
        __m128 t3 = _mm_mul_ps(v[3], vec.v);
        _MM_TRANSPOSE4_PS(t0, t1, t2, t3);
        return vec3sse(_mm_add_ps(_mm_add_ps(t0, t1), _mm_add_ps(t2, t3)));
    }

    matrix4sse operator*(matrix4sse rhs) const {
        _MM_TRANSPOSE4_PS(rhs.v[0], rhs.v[1], rhs.v[2], rhs.v[3]);
        __m128 r[4];
        for (int i = 0; i < 4; ++i) {
            __m128 t0 = _mm_mul_ps(v[0], rhs.v[i]);
            __m128 t1 = _mm_mul_ps(v[1], rhs.v[i]);
            __m128 t2 = _mm_mul_ps(v[2], rhs.v[i]);
            __m128 t3 = _mm_mul_ps(v[3], rhs.v[i]);
            _MM_TRANSPOSE4_PS(t0, t1, t2, t3);
            r[i] = _mm_add_ps(_mm_add_ps(t0, t1), _mm_add_ps(t2, t3));
        }
        _MM_TRANSPOSE4_PS(r[0], r[1], r[2], r[3]); //fixme
        return matrix4sse(r[0], r[1], r[2], r[3]);
    }

    __m128 v[4];
};

}

#endif  // OSDUTIL_MATH_SSE_H
