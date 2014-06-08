#ifndef OSDUTIL_BEZIER_H
#define OSDUTIL_BEZIER_H

//#include "../../../opensubdiv/version.h"  // FIXME
#include <cmath>

//namespace OpenSubdiv {
//namespace OPENSUBDIV_VERSION {

struct float3 {
    typedef float ElementType;
    static const int LENGTH = 3;

    float3(float f) {
        v[0] = v[1] = v[2] = f;
    }
    float3() {
        v[0] = v[1] = v[2] = 0;
    }
    float3(float x, float y, float z) {
        v[0] = x; v[1] = y; v[2] = z;
    }
    float3 operator*(float f) const { return float3(v[0] * f, v[1] * f, v[2] * f); }
    float3 operator+(float3 f) const { return float3(v[0] + f.v[0], v[1] + f.v[1], v[2] + f.v[2]); }
    float3 operator-(float3 f) const { return float3(v[0] - f.v[0], v[1] - f.v[1], v[2] - f.v[2]); }
    float3 &operator+=(const float3 &f2) { v[0] += f2.v[0]; v[1] += f2.v[1]; v[2] += f2.v[2]; return (*this); }
    float3 &operator-=(const float3 &f2) { v[0] -= f2.v[0]; v[1] -= f2.v[1]; v[2] -= f2.v[2]; return (*this); }

    float operator[](int i) const { return v[i]; }
    float &operator[](int i) { return v[i]; }
    float v[3];

    float length() const { return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }
    void normalize() {
        float len = length();
        //    if (fabs(len) > 1.0e-6) {
        float inv_len = 1.0 / len;
        v[0] *= inv_len;
        v[1] *= inv_len;
        v[2] *= inv_len;
        //    }
    }
};

inline static float3 operator*(float f, const float3 &v) {
    return v*f;
}
inline float3 cross(float3 a, float3 b) {
    float3 c;
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
    return c;
}

// ----------------------------------------------------------------------------

struct Matrix4 {
    Matrix4() { }
    Matrix4(float m00, float m01, float m02, float m03,
            float m10, float m11, float m12, float m13,
            float m20, float m21, float m22, float m23,
            float m30, float m31, float m32, float m33) {
            m[0] = m00;
            m[1] = m01;
            m[2] = m02;
            m[3] = m03;

            m[4] = m10;
            m[5] = m11;
            m[6] = m12;
            m[7] = m13;

            m[8] = m20;
            m[9] = m21;
            m[10] = m22;
            m[11] = m23;

            m[12] = m30;
            m[13] = m31;
            m[14] = m32;
            m[15] = m33;
    }

    float * operator[] (int i){            //pointer what self is const.
        return element[i];
    }
    const float * operator[] (int i)const{//pointer what self is const.
        return element[i];
    }

    union {
        float m[16];
        float element[4][4];
    };
};

static    inline Matrix4 operator*(const Matrix4& lhs, const Matrix4& rhs)  {
    	return Matrix4(
            lhs[0][0] * rhs[0][0] + lhs[0][1] * rhs[1][0] + lhs[0][2] * rhs[2][0] + lhs[0][3] * rhs[3][0],
            lhs[0][0] * rhs[0][1] + lhs[0][1] * rhs[1][1] + lhs[0][2] * rhs[2][1] + lhs[0][3] * rhs[3][1],
            lhs[0][0] * rhs[0][2] + lhs[0][1] * rhs[1][2] + lhs[0][2] * rhs[2][2] + lhs[0][3] * rhs[3][2],
            lhs[0][0] * rhs[0][3] + lhs[0][1] * rhs[1][3] + lhs[0][2] * rhs[2][3] + lhs[0][3] * rhs[3][3],
        
            lhs[1][0] * rhs[0][0] + lhs[1][1] * rhs[1][0] + lhs[1][2] * rhs[2][0] + lhs[1][3] * rhs[3][0],
            lhs[1][0] * rhs[0][1] + lhs[1][1] * rhs[1][1] + lhs[1][2] * rhs[2][1] + lhs[1][3] * rhs[3][1],
            lhs[1][0] * rhs[0][2] + lhs[1][1] * rhs[1][2] + lhs[1][2] * rhs[2][2] + lhs[1][3] * rhs[3][2],
            lhs[1][0] * rhs[0][3] + lhs[1][1] * rhs[1][3] + lhs[1][2] * rhs[2][3] + lhs[1][3] * rhs[3][3],
        
            lhs[2][0] * rhs[0][0] + lhs[2][1] * rhs[1][0] + lhs[2][2] * rhs[2][0] + lhs[2][3] * rhs[3][0],
            lhs[2][0] * rhs[0][1] + lhs[2][1] * rhs[1][1] + lhs[2][2] * rhs[2][1] + lhs[2][3] * rhs[3][1],
            lhs[2][0] * rhs[0][2] + lhs[2][1] * rhs[1][2] + lhs[2][2] * rhs[2][2] + lhs[2][3] * rhs[3][2],
            lhs[2][0] * rhs[0][3] + lhs[2][1] * rhs[1][3] + lhs[2][2] * rhs[2][3] + lhs[2][3] * rhs[3][3],
        
            lhs[3][0] * rhs[0][0] + lhs[3][1] * rhs[1][0] + lhs[3][2] * rhs[2][0] + lhs[3][3] * rhs[3][0],
            lhs[3][0] * rhs[0][1] + lhs[3][1] * rhs[1][1] + lhs[3][2] * rhs[2][1] + lhs[3][3] * rhs[3][1],
            lhs[3][0] * rhs[0][2] + lhs[3][1] * rhs[1][2] + lhs[3][2] * rhs[2][2] + lhs[3][3] * rhs[3][2],
            lhs[3][0] * rhs[0][3] + lhs[3][1] * rhs[1][3] + lhs[3][2] * rhs[2][3] + lhs[3][3] * rhs[3][3]        
        );
    }

    inline float3 operator*(const Matrix4& m, const float3& v)
    {
        float r[4];
        r[0] = m.m[0] * v[0] + m.m[1] * v[1] + m.m[2] * v[2] + m.m[3];
        r[1] = m.m[4] * v[0] + m.m[5] * v[1] + m.m[6] * v[2] + m.m[7];
        r[2] = m.m[8] * v[0] + m.m[9] * v[1] + m.m[10] * v[2] + m.m[11];
        r[3] = m.m[12] * v[0] + m.m[13] * v[1] + m.m[14] * v[2] + m.m[15];
        
        float ir = 1.0f/r[3];
        return float3(r[0]*ir,r[1]*ir,r[2]*ir);
    }
// ----------------------------------------------------------------------------

template<class T>
class matrix3t {
public:
    matrix3t() { }
    matrix3t(T m00, T m01, T m02, T m10, T m11, T m12, T m20, T m21, T m22) {
        m[0] = m00;
        m[1] = m01;
        m[2] = m02;
        m[3] = m10;
        m[4] = m11;
        m[5] = m12;
        m[6] = m20;
        m[7] = m21;
        m[8] = m22;
    }
    matrix3t(const matrix3t& rhs){
        memcpy(m,rhs.m,sizeof(T)*9);
    }
    T* operator[](size_t i)     {return m+3*i;}
    const T* operator[](size_t i)const{return m+3*i;}

    //0,1,2
    //3,4,5
    //6.7.8
    void transpose(){
        std::swap(m[1],m[3]);
        std::swap(m[2],m[6]);
        std::swap(m[5],m[7]);
    }
    T m[3*3];
};
typedef matrix3t<float> matrix3x;

template<class T, class V>
inline V operator*(const matrix3t<T>& m, const V& v)
{
    V r;
    r[0] = m.m[0] * v[0] + m.m[1] * v[1] + m.m[2] * v[2];
    r[1] = m.m[3] * v[0] + m.m[4] * v[1] + m.m[5] * v[2];
    r[2] = m.m[6] * v[0] + m.m[7] * v[1] + m.m[8] * v[2];
    return r;
}

// ----------------------------------------------------------------------------

template <class T, int N>
inline void transpose_matrix_n(T m[])
{
    for (int y = 0; y < N; ++y) {
        for (int x = y+1; x < N; ++x) {
            std::swap(m[x*N+y], m[y*N+x]);
        }
    }
}

template<int N, int k>
struct Binomial
{
    const static int value =  (Binomial<N-1,k-1>::value + Binomial<N-1,k>::value);
};

template<>
struct Binomial<0, 0>
{
    const static int value = 1;
};

template<int N>
struct Binomial<N, 0>
{
    const static int value = 1;
};

template<int N>
struct Binomial<N, N>
{
    const static int value = 1;
};

template<int N, int k>
inline int binomial()
{
    return Binomial<N, k>::value;
}

template<typename REAL, int N, int I>
inline REAL bernstein(REAL t)
{
    return REAL(binomial<N, I>()) * std::pow(REAL(1)-t, N-I) * std::pow(t, I);
}

template<typename T, typename REAL, int n, int i>
struct BezierCurveEvaluate
{
   static T evaluate(const REAL t, const T *const cp)
   {
       return cp[i] * bernstein<REAL, n, i>(t) +
          BezierCurveEvaluate<T, REAL, n, i - 1>::evaluate(t, cp);
   }
};

template<typename T, typename REAL, int N>
struct BezierCurveEvaluate<T, REAL, N, 0>
{
   static T evaluate(const REAL t, const T *const cp)
   {
       return cp[0] * bernstein<REAL, N, 0>(t);
   }
};

template<typename T, typename REAL, int N>
T evaluate(REAL t, const T * const cp)
{
    return BezierCurveEvaluate<T, REAL, N-1, N-1>::evaluate(t, cp);
}
template<typename T, typename REAL, int N>
T evaluateD(REAL t, const T * const cp)
{
    if (N == 4) {
        REAL t2 = t*t;
        return cp[0] * (3*t2*-1 + 2*t* 3 + -3)
            + cp[1] * (3*t2* 3 + 2*t*-6 +  3)
            + cp[2] * (3*t2*-3 + 2*t* 3)
            + cp[3] * (3*t2* 1);
    }
}

template <int N, int I, typename VALUE_TYPE, typename REAL>
struct BezierSplit {
    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    BezierSplit(ValueType a[], ValueType b[], const ValueType p[], Real t) {
        int tn = N-1-I;//3
        a[I] = p[0];
        b[tn] = p[tn];
        ValueType tmp[N-1];
        for (int j = 0; j < tn; j++) {
            tmp[j] = (1-t)*p[j] + t*p[j+1];//p[j] + t*(p[j+1] - p[j])
        }
        BezierSplit<N, I+1, VALUE_TYPE, REAL>(a, b, &tmp[0], t);
    }
};
template <int N, typename VALUE_TYPE, typename REAL>
struct BezierSplit<N, N, VALUE_TYPE, REAL> {
    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    BezierSplit(ValueType[], ValueType[], const ValueType[], Real) {}
};

template<int N=4, class VALUE_TYPE=float3, class REAL=float>
class OsdUtilBezierPatch {
public:
    typedef OsdUtilBezierPatch<N, VALUE_TYPE, REAL> This;
    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    static const int Order = N;
    static const int Ncp = N*N;

    /// Constructor.
    OsdUtilBezierPatch() {
        //for (int i = 0; i < ncp; ++i) _cp[i] = VALUE_TYPE();
    }

    /// Constructor.
    OsdUtilBezierPatch(OsdUtilBezierPatch const &other) {
        for (int i = 0; i < Ncp; ++i) _cp[i] = other._cp[i];
    }

    /// Constructor.
    template<class T>
    OsdUtilBezierPatch(T const *p) {
        for (int i = 0; i < Ncp; ++i) _cp[i] = p[i];
    }

    /// Get
    ValueType const &Get(int i, int j) const {
        return _cp[N*j + i];
    }

    /// Set
    void Set(int i, int j, VALUE_TYPE const &v) {
        _cp[N*j + i] = v;
    }

    /// Evaluate
    ValueType Evaluate(REAL u, REAL v) const {
        ValueType b[N];
        for (int i = 0; i < N; ++i) {
            b[i] = evaluate<ValueType, Real, N>(u, _cp+i*N);
        }
        return evaluate<ValueType, Real, N>(v, b);
    }

    /// Evaluate derivative
    ValueType EvaluateDu(REAL u, REAL v) const {
        ValueType b[N];
        for (int i = 0; i < N; ++i) {
            b[i] = evaluateD<ValueType, Real, N>(u, _cp+i*N);
        }
        return evaluate<ValueType, Real, N>(v, b);
    }

    /// Evaluate derivative
    ValueType EvaluateDv(REAL u, REAL v) const {
        ValueType b[N];
        for (int i = 0; i < N; ++i) {
            b[i] = evaluate<ValueType, Real, N>(u, _cp+i*N);
        }
        return evaluateD<ValueType, Real, N>(v, b);
    }

    /// min
    ValueType GetMin() const {
        ValueType min = _cp[0];
        for (int i = 1; i < Ncp; ++i) {
            for (int j = 0; j < ValueType::LENGTH; ++j) {
                min[j] = std::min(min[j], _cp[i][j]);
            }
        }
        return min;
    }
    /// max
    ValueType GetMax() const {
        ValueType max = _cp[0];
        for (int i = 1; i < Ncp; ++i) {
            for (int j = 0; j < ValueType::LENGTH; ++j) {
                max[j] = std::max(max[j], _cp[i][j]);
            }
        }
        return max;
    }
    void GetMinMax(ValueType &min, ValueType &max) const {
        min = max = _cp[0];
        for (int i = 1; i < Ncp; ++i) {
            for (int j = 0; j < ValueType::LENGTH; ++j) {
                min[j] = std::min(min[j], _cp[i][j]);
                max[j] = std::max(max[j], _cp[i][j]);
            }
        }
    }
    void GetMinMax(ValueType &min, ValueType &max, Real eps) const {
        GetMinMax(min, max);
        min -= ValueType(eps);
        max += ValueType(eps);
    }

    /// transform
    template<typename MATRIX>
    This & Transform(const MATRIX &m) {
        for (int i = 0; i < Ncp; ++i) _cp[i] = m * _cp[i];
        return *this;
    }

    /// split
    void Split(This patches[4], Real u, Real v) const {
        ValueType tmp0[Ncp], tmp1[Ncp];

        // split U
        for (int i = 0; i < N; ++i) {
            bezierSplit(&tmp0[i*N+0], &tmp1[i*N+0], &_cp[i*N+0], u);
        }
        // uv -> vu
        transpose_matrix_n<ValueType, N>(&tmp0[0]);// 00 01
        transpose_matrix_n<ValueType, N>(&tmp1[0]);// 10 11

        // split V
        for (int i = 0; i < N; ++i) {
            bezierSplit(&patches[0]._cp[i*N+0],
                        &patches[2]._cp[i*N+0],
                        &tmp0[i*N+0], v);
            bezierSplit(&patches[1]._cp[i*N+0],
                        &patches[3]._cp[i*N+0],
                        &tmp1[i*N+0], v);
        }
        // vu -> uv
        transpose_matrix_n<ValueType, N>(patches[0]._cp);//00
        transpose_matrix_n<ValueType, N>(patches[1]._cp);//10
        transpose_matrix_n<ValueType, N>(patches[2]._cp);//01
        transpose_matrix_n<ValueType, N>(patches[3]._cp);//11
    }

    /// split
    void SplitU(This patches[2], Real u) const {
        for (int i = 0; i < N; ++i) {
            bezierSplit(&patches[0]._cp[i*N+0],
                        &patches[1]._cp[i*N+0],
                        &_cp[i*N+0], u);
        }
    }

    /// split
    void SplitV(This patches[2], Real v) const {
        This tmp(*this);
        transpose_matrix_n<ValueType, N>(tmp._cp);
        for (int i = 0; i < N; ++i) {
            bezierSplit(&patches[0]._cp[i*N+0],
                        &patches[1]._cp[i*N+0],
                        &_cp[i*N+0], v);
        }
        transpose_matrix_n<ValueType, N>(patches[0]._cp);
        transpose_matrix_n<ValueType, N>(patches[1]._cp);
    }

    /// crop
    void CropU(This &patch, Real u0, Real u1) const {
        for (int i = 0; i < N; ++i) {
            ValueType tmp0[N], tmp1[N];
            bezierSplit(&tmp0[0], &tmp1[0], &_cp[i*N+0], u1);
            bezierSplit(&tmp1[0], &patch._cp[i*N+0], &tmp0[0],  u0/u1);
        }
    }

    /// crop
    void CropV(This &patch, Real v0, Real v1) const {
        This tmp(*this);
        transpose_matrix_n<ValueType, N>(tmp._cp);
        for (int i = 0; i < N; ++i) {
            ValueType tmp0[N], tmp1[N];
            bezierSplit(&tmp0[0], &tmp1[0], &tmp._cp[i*N+0], v1);
            bezierSplit(&tmp1[0], &patch._cp[i*N+0], &tmp0[0],  v0/v1);
        }
        transpose_matrix_n<ValueType, N>(patch._cp);
    }

    void Test() {
        This patches[4];
        Split(patches, 0.5, 0.5);
        SplitU(patches, 0.5);
        SplitV(patches, 0.5);
        CropU(patches[0], 0.25, 0.75);
        CropV(patches[0], 0.25, 0.75);
    }

    ValueType GetLv() const {
        return Get(0, N-1) - Get(0, 0) + Get(N-1, N-1) - Get(N-1, 0);
    }

    ValueType GetLu() const {
        return Get(N-1, 0) - Get(0, 0) + Get(N-1, N-1) - Get(0, N-1);
    }

private:
    void bezierSplit(ValueType a[], ValueType b[], const ValueType p[], Real t) const {
        BezierSplit<N, 0, ValueType, Real> split(a, b, p, t);
    }

    ValueType _cp[Ncp];
};


//}  // end namespace OPENSUBDIV_VERSION
//using namespace OPENSUBDIV_VERSION;

//}  // end namespace OpenSubdiv

#endif  // OSDUTIL_BEZIER_H
