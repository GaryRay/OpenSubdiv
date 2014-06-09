#ifndef OSDUTIL_BEZIER_H
#define OSDUTIL_BEZIER_H

#include <cmath>
#include "math.h"

//#include "../../../opensubdiv/version.h"  // FIXME
//namespace OpenSubdiv {
//namespace OPENSUBDIV_VERSION {

namespace OsdUtil {

// ----------------------------------------------------------------------------
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
inline REAL bernstein(REAL t) {
    return REAL(binomial<N, I>() * std::pow(REAL(1)-t, N-I) * std::pow(t, I));
}

template<typename T, typename REAL, int N, int I>
struct BezierCurveEvaluate
{
    static T evaluate(const REAL t, const T *const cp) {
        return cp[I] * bernstein<REAL, N, I>(t) +
            BezierCurveEvaluate<T, REAL, N, I - 1>::evaluate(t, cp);
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

template <int N, int I, int STRIDE, typename VALUE_TYPE, typename REAL>
struct BezierSplit {
    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    BezierSplit(ValueType a[], ValueType b[], const ValueType p[], Real t) {
        int tn = N-1-I;//3
        a[I*STRIDE] = p[0];
        b[tn*STRIDE] = p[tn*STRIDE];
        ValueType tmp[(N-1)*STRIDE];
        for (int j = 0; j < tn; j++) {
            tmp[j*STRIDE] = (1-t)*p[j*STRIDE] + t*p[(j+1)*STRIDE];//p[j] + t*(p[j+1] - p[j])
        }
        BezierSplit<N, I+1, STRIDE, VALUE_TYPE, REAL>(a, b, &tmp[0], t);
    }
};

template <int N, int STRIDE, typename VALUE_TYPE, typename REAL>
struct BezierSplit<N, N, STRIDE, VALUE_TYPE, REAL> {
    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    BezierSplit(ValueType[], ValueType[], const ValueType[], Real) {
    }
};

template<class VALUE_TYPE, class REAL, int N=4>
class OsdUtilBezierPatch {
public:
    typedef OsdUtilBezierPatch<VALUE_TYPE, REAL, N> This;
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
            b[i] = evaluate(u, _cp+i*N);
        }
        return evaluate(v, b);
    }

    /// Evaluate derivative
    ValueType EvaluateDu(REAL u, REAL v) const {
        ValueType b[N];
        for (int i = 0; i < N; ++i) {
            b[i] = evaluateD(u, _cp+i*N);
        }
        return evaluate(v, b);
    }

    /// Evaluate derivative
    ValueType EvaluateDv(REAL u, REAL v) const {
        ValueType b[N];
        for (int i = 0; i < N; ++i) {
            b[i] = evaluate(u, _cp+i*N);
        }
        return evaluateD(v, b);
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
        transpose(&tmp0[0]);// 00 01
        transpose(&tmp1[0]);// 10 11

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
        transpose(patches[0]._cp);//00
        transpose(patches[1]._cp);//10
        transpose(patches[2]._cp);//01
        transpose(patches[3]._cp);//11
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
        for (int i = 0; i < N; ++i) {
            bezierSplitV(&patches[0]._cp[i],
                         &patches[1]._cp[i],
                         &_cp[i], v);
        }
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
        for (int i = 0; i < N; ++i) {
            ValueType tmp0[N*N], tmp1[N*N];
            bezierSplitV(&tmp0[i], &tmp1[i], &_cp[i], v1);
            bezierSplitV(&tmp1[i], &patch._cp[i], &tmp0[i],  v0/v1);
        }
    }

    ValueType GetLv() const {
        return Get(0, N-1) - Get(0, 0) + Get(N-1, N-1) - Get(N-1, 0);
    }

    ValueType GetLu() const {
        return Get(N-1, 0) - Get(0, 0) + Get(N-1, N-1) - Get(0, N-1);
    }

    void Test() {
        This patches[4];
        Split(patches, 0.5, 0.5);
        SplitU(patches, 0.5);
        SplitV(patches, 0.5);
        CropU(patches[0], 0.25, 0.75);
        CropV(patches[0], 0.25, 0.75);
    }

private:
    void bezierSplit(ValueType a[], ValueType b[], const ValueType p[], Real t) const {
        BezierSplit<N, 0, /*stride*/1, ValueType, Real> split(a, b, p, t);
    }
    void bezierSplitV(ValueType a[], ValueType b[], const ValueType p[], Real t) const {
        BezierSplit<N, 0, /*stride*/N, ValueType, Real> split(a, b, p, t);
    }

    static void transpose(ValueType m[]) {
        for (int y = 0; y < N; ++y) {
            for (int x = y+1; x < N; ++x) {
                std::swap(m[x*N+y], m[y*N+x]);
            }
        }
    }

    static ValueType evaluate(Real t, const ValueType * const cp) {
        return BezierCurveEvaluate<ValueType, Real, N-1, N-1>::evaluate(t, cp);
    }
    static ValueType evaluateD(Real t, const ValueType * const cp) {
        if (N == 4) {
            REAL t2 = t*t;
            return cp[0] * (3*t2*-1 + 2*t* 3 + -3)
                + cp[1] * (3*t2* 3 + 2*t*-6 +  3)
                + cp[2] * (3*t2*-3 + 2*t* 3)
                + cp[3] * (3*t2* 1);
        } else {
            return ValueType(0);
            // not implemented !
        }
    }

    ValueType _cp[Ncp];
};

}   // end namespace OsdUtil

//}  // end namespace OPENSUBDIV_VERSION
//using namespace OPENSUBDIV_VERSION;

//}  // end namespace OpenSubdiv

#endif  // OSDUTIL_BEZIER_H
