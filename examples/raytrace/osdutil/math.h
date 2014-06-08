#ifndef OSDUTIL_MATH_H
#define OSDUTIL_MATH_H

//#include "../../../opensubdiv/version.h"  // FIXME
//namespace OpenSubdiv {
//namespace OPENSUBDIV_VERSION {

namespace OsdUtil {

template <typename REAL>
struct vec3 {
    typedef REAL Real;
    typedef REAL ElementType;
    static const int LENGTH = 3;

    vec3(Real f) {
        v[0] = v[1] = v[2] = f;
    }
    vec3() {
        v[0] = v[1] = v[2] = 0;
    }
    vec3(Real x, Real y, Real z) {
        v[0] = x; v[1] = y; v[2] = z;
    }
    vec3 operator*(Real f) const {
        return vec3(v[0] * f, v[1] * f, v[2] * f);
    }
    vec3 operator+(vec3 f) const {
        return vec3(v[0] + f.v[0], v[1] + f.v[1], v[2] + f.v[2]);
    }
    vec3 operator-(vec3 f) const {
        return vec3(v[0] - f.v[0], v[1] - f.v[1], v[2] - f.v[2]);
    }
    vec3 &operator+=(const vec3 &f2) {
        v[0] += f2.v[0]; v[1] += f2.v[1]; v[2] += f2.v[2];
        return (*this);
    }
    vec3 &operator-=(const vec3 &f2) {
        v[0] -= f2.v[0]; v[1] -= f2.v[1]; v[2] -= f2.v[2];
        return (*this);
    }
    Real operator[](int i) const { return v[i]; }
    Real &operator[](int i) { return v[i]; }

    Real length() const { return std::sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }
    void normalize() {
        Real inv_len = 1.0 / length();
        v[0] *= inv_len;
        v[1] *= inv_len;
        v[2] *= inv_len;
    }

    float v[3];
};

template <typename REAL>
inline static vec3<REAL> operator*(REAL f, const vec3<REAL> &v) {
    return v*f;
}

template <typename REAL>
inline static vec3<REAL> cross(vec3<REAL> const &a, vec3<REAL> const &b) {
    vec3<REAL> c;
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

    template <typename REAL>
    inline vec3<REAL> operator*(const Matrix4& m, vec3<REAL> const & v)
    {
        REAL r[4];
        r[0] = m.m[0] * v[0] + m.m[1] * v[1] + m.m[2] * v[2] + m.m[3];
        r[1] = m.m[4] * v[0] + m.m[5] * v[1] + m.m[6] * v[2] + m.m[7];
        r[2] = m.m[8] * v[0] + m.m[9] * v[1] + m.m[10] * v[2] + m.m[11];
        r[3] = m.m[12] * v[0] + m.m[13] * v[1] + m.m[14] * v[2] + m.m[15];
        REAL ir = REAL(1.0)/r[3];
        return vec3<REAL>(r[0]*ir, r[1]*ir, r[2]*ir);
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

template<class T, class V>
inline V operator*(const matrix3t<T>& m, const V& v)
{
    V r;
    r[0] = m.m[0] * v[0] + m.m[1] * v[1] + m.m[2] * v[2];
    r[1] = m.m[3] * v[0] + m.m[4] * v[1] + m.m[5] * v[2];
    r[2] = m.m[6] * v[0] + m.m[7] * v[1] + m.m[8] * v[2];
    return r;
}

//----------------------------------------------------------

typedef matrix3t<float> matrix3x;
typedef vec3<float> float3;

}

#endif  // OSDUTIL_MATH_H
