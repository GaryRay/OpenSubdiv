#ifndef OSDUTIL_BILINEAR_INTERSECT_H
#define OSDUTIL_BILINEAR_INTERSECT_H

#include <limits>

template<class T>
static bool solveBilinearPatch(T* t, T* u, T* v, T tmin, T tmax, T tt, T uu, T vv) {
    if (tmin <= tt && tt <= tmax) {
        *t = tt;
        *u = uu;
        *v = vv;
        return true;
    }
    return false;
}

template <class T>
static int solve2(T root[2], T const coeff[3]) {
    T A = coeff[0];
    T B = coeff[1];
    T C = coeff[2];
    if (fabs(A) <= std::numeric_limits<T>::min()) {
        if (fabs(B) <= std::numeric_limits<T>::min()) return 0;
        T x = -C/B;
        root[0] = x;
        return 1;
    } else {
        T D = B*B - 4*A*C;
        if (D < 0) {
            return 0;
        } else if (D == 0) {
            T x = -0.5*B/A;
            root[0] = x;
            return 1;
        } else {
            T x1 = (fabs(B) + sqrt(D)) / (2.0 * A);
            if (B >= 0) {
                x1 = -x1;
            }
            T x2 = C / (A * x1);
            if (x1 > x2) std::swap(x1, x2);
            root[0] = x1;
            root[1] = x2;
            //            printf ("(%f, %f), ", x1, x2);
            return 2;
        }
    }
}
template<class T>
static T computeU(T A1, T A2, T B1, T B2, T C1, T C2, T D1, T D2, T v) {
    //return div((v*(C1-C2)+(D1-D2)),(v*(A2-A1)+(B2-B1)));
    T a = v*A2 + B2;
    T b = v*(A2-A1) + B2 - B1;
    if (fabs(b) >= fabs(a)) {
        return (v*(C1-C2) + D1 - D2)/b;
    } else {
        return (-v*C2 - D2)/a;
    }
}
template<class T>
static T computeT(T a, T b, T c, T d, T iq, T u, T v)
{
    return ((u*v)*a + u*b + v*c + d)*iq;
}

template<typename T, typename V>
static bool testBilinearPatch(T *t, T *u, T *v,
                              V const p[4], T tmin, T tmax) {
    typedef typename V::ElementType Scalar;

    V const & p00 = p[0];
    V const & p10 = p[1];
    V const & p01 = p[2];
    V const & p11 = p[3];

    static const int nPlane = 2;
    V a = p11-p10-p01+p00;
    V b = p10-p00;
    V c = p01-p00;
    V d = p00;

    //xz-zx
    Scalar A1 = a[0];
    Scalar B1 = b[0];
    Scalar C1 = c[0];
    Scalar D1 = d[0];

    //yz-zy
    Scalar A2 = a[1];
    Scalar B2 = b[1];
    Scalar C2 = c[1];
    Scalar D2 = d[1];

    Scalar F1 = A2*C1-A1*C2;
    Scalar F2 = A2*D1-A1*D2+B2*C1-B1*C2;
    Scalar F3 = B2*D1-B1*D2;

    Scalar coeff[] = {F1,F2,F3};
    Scalar root[2] = {};
    //    printf("<%f, %f, %f>\n", F1, F2, F3);
    int nRet = solve2(root,coeff);
    if (nRet) {
        bool bRet = false;
        for (int i = 0; i < nRet; ++i) {
            Scalar vv = root[i];
            if (0 <= vv && vv <= 1) {
                Scalar uu = computeU(A1, A2, B1, B2, C1, C2, D1, D2, vv);
                if (0 <= uu && uu <= 1) {
                    Scalar tt = computeT(a[nPlane], b[nPlane], c[nPlane], d[nPlane],
                                         Scalar(1), uu, vv);
                    if (solveBilinearPatch(t, u, v, tmin, tmax, tt, uu, vv)){
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


#endif  // OSDUTIL_BILINEAR_INTERSECT_H
