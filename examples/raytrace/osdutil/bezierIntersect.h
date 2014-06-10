#ifndef OSDUTIL_BEZIER_INTERSECT_H
#define OSDUTIL_BEZIER_INTERSECT_H

#include "bezier.h"
#include "bilinearIntersect.h"

#include <limits>
#include <cfloat>
#include "../common.h"

//namespace OpenSubdiv {
//namespace OPENSUBDIV_VERSION {

#define DIRECT_BILINEAR 1
#define USE_BEZIERCLIP 1
#define USE_COARSESORT 1

#ifdef __GNUC__
#define NO_INLINE __attribute__((noinline))
#else
#define NO_INLINE
#endif

namespace OsdUtil {

template <typename REAL>
struct Epsilon {
    static const REAL EPS = 1e-4;
    static const REAL EPS2 = FLT_MIN;
};
template <>
struct Epsilon<double> {
    static const double EPS = 1e-16;    // revisit
    static const double EPS2 = 1e-37;
};

template<class VALUE_TYPE, class REAL, int N=4>
class OsdUtilBezierPatchIntersection {
public:

    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    typedef OsdUtilBezierPatch<ValueType, Real, N> PatchType;

    static const REAL EPS = Epsilon<Real>::EPS;
    static const REAL EPS2 = Epsilon<Real>::EPS2;
    static const REAL UVEPS = 1.0/32.0;
    static const int  DEFAULT_MAX_LEVEL = 10;

    struct RangeAABB {
        Real tmin, tmax;
    };
    struct UVT {
        Real u, v, t;
        int level;
    };

    OsdUtilBezierPatchIntersection(PatchType const &patch, int maxLevel=DEFAULT_MAX_LEVEL) NO_INLINE :
        _patch(patch), _maxLevel(maxLevel) {
        _uRange[0] = _vRange[0] = 0;
        _uRange[1] = _vRange[1] = 1;

        patch.GetMinMax(_min, _max, 0.01);
        _eps = EPS;
    }
    ~OsdUtilBezierPatchIntersection() { }

    bool Test(Intersection* info, Ray const &r, Real tmin, Real tmax) const NO_INLINE {

        RangeAABB rng;
        if (intersectAABB(&rng, _min, _max, r, tmin, tmax)) {
            tmin = std::max(tmin, rng.tmin);
            tmax = std::min(tmax, rng.tmax);

            return testInternal(info, r, tmin, tmax);
        }
        return false;
    }

protected:
    template <typename MATRIX4>
    void getZAlign(MATRIX4 & mat, const Ray &r) const {
        ValueType org(r.org[0], r.org[1], r.org[2]);
        ValueType dir(r.dir[0], r.dir[1], r.dir[2]);
        ValueType z = dir;

        int plane = 0;
        if (fabs(z[1]) < fabs(z[plane])) plane = 1;
        if (fabs(z[2]) < fabs(z[plane])) plane = 2;

        ValueType x = (plane == 0) ? ValueType(1, 0, 0) :
            ((plane == 1) ? ValueType(0, 1, 0) : ValueType(0, 0, 1));
        ValueType y = cross(z,x);
        y.normalize();
        x = cross(y,z);
        MATRIX4 rot = MATRIX4(x[0],x[1],x[2],0,
                              y[0],y[1],y[2],0,
                              z[0],z[1],z[2],0,
                              0,0,0,1);
        MATRIX4 trs = MATRIX4(1,0,0,-org[0],
                              0,1,0,-org[1],
                              0,0,1,-org[2],
                              0,0,0,1);
        mat = rot*trs;
    }

    bool testInternal(Intersection* info, const Ray& r, Real tmin, Real tmax) const NO_INLINE {
        typename ValueType::Matrix4Type mat;
        getZAlign(mat, r);

        UVT uvt;
        PatchType patch(_patch);
        patch.Transform(mat);
        if (testBezierPatch(&uvt, patch, tmin, tmax, _eps)) {
            Real t = uvt.t;
            Real u = uvt.u;
            Real v = uvt.v;

            u = _uRange[0]*(1-u) + _uRange[1]*u;//global
            v = _vRange[0]*(1-v) + _vRange[1]*v;//global
            info->t = t;
            info->u = u;
            info->v = v;
            info->clipLevel = uvt.level;
            {
                ValueType du = _patch.EvaluateDu(u,v);
                ValueType dv = _patch.EvaluateDv(u,v);
                ValueType normal = cross(du,dv);
                normal.normalize();
                info->normal = real3(normal[0], normal[1], normal[2]);
                info->geometricNormal = real3(normal[0], normal[1], normal[2]);
                //                info->tangent  = Conv(U);
                //                info->binormal = Conv(V);
            }
            return true;
        }
        return false;
    }
    bool testBezierPatch(UVT* info, PatchType const & patch, Real zmin, Real zmax, Real eps) const {
        ValueType min, max;
        patch.GetMinMax(min, max, EPS*1e-3);

        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2]) return false;//z

        return testBezierClipU(info, patch, 0, 1, 0, 1, zmin, zmax, 0, _maxLevel, eps);
    }

    // ----------------------------------------------------------------------
    template <typename MATRIX3>
    static void getRotate(MATRIX3 & mat, ValueType const &dx) {
        const ValueType z(0, 0, 1);
        ValueType x = dx;
        typename ValueType::ElementType l = 1.0/sqrt(x[0]*x[0] + x[1]*x[1]);
        x = ValueType(x[0]*l, x[1]*l, 0);
        ValueType y(-x[1], x[0], 0);
        mat = MATRIX3(x[0],x[1],x[2],
                      y[0],y[1],y[2],
                      z[0],z[1],z[2]);
    }

    static void rotateU(PatchType& patch)
    {
        ValueType dx = patch.GetLv();
        typename ValueType::Matrix3Type rot;
        getRotate(rot, dx);
        patch.Transform(rot);
    }

    static void rotateV(PatchType& patch)
    {
        ValueType dx = patch.GetLu();
        typename ValueType::Matrix3Type rot;
        getRotate(rot, dx);
        patch.Transform(rot);
    }

    static bool isEps(ValueType const & min, ValueType const & max, Real eps) {
        //float zw = max[2]-min[2];
        //if(zw<=eps)return true;
        //REAL xd = std::max<REAL>(fabs(min[0]),max[0]);
        //REAL yd = std::max<REAL>(fabs(min[1]),max[1]);
        //if(!(xd<=eps&&yd<=eps))return false;

        REAL xw = max[0]-min[0];
        REAL yw = max[1]-min[1];
        //
        if(xw<=eps && yw<=eps)return true;
        else return false;
    }

    static bool isLevel(int level, int max_level) {
        return (level>=max_level);
        //return false;
    }

    static bool isClip(int level) {
        static const int div_level[] =
        {
            1,
            1,1,
            2,2,//4
            3,3,3,3,//8
            4,4,4,4,4,4,4,4,//16
        };
        int nlevel = 0;
        if (N<=16) {
            nlevel = div_level[N];
        }else{
            nlevel = (int)ceil(log2(N));
        }
        return nlevel*2<=level;
    }
    static Real lerp(Real a, Real b, Real t) {
        return a+(b-a)*t;
    }

    static void coarseSort(int order[2], PatchType tmp[2])
    {
#if USE_COARSESORT
        Real zs[2];
        for (int i = 0; i < 2; ++i) {
            int u = PatchType::Order>>1;
            int v = PatchType::Order>>1;

            zs[i] = tmp[i].Get(u,v)[2];
        }
        if (zs[0] <= zs[1]) {
            order[0] = 0;
            order[1] = 1;
        } else {
            order[0] = 1;
            order[1] = 0;
        }
#endif
    }

    // --------------------------------------------------------------------

    template<typename VECTOR>
    static bool scanMinMax(VECTOR const & p) {
        int nUpper = 0, nLower = 0;
        for (int i = 0; i < VECTOR::LENGTH; ++i) {
            if (p[i][1] > 0) nUpper++;
            else nLower++;
            if (nUpper && nLower) return true;
        }
        return false;
    }

    template<typename T>
    static T slope(T const a[], T const b[]) {
        T d0 = b[0] - a[0];
        T d1 = b[1] - a[2];
        return fabs(d0 / d1);
    }

    template<typename T>
    static T dist(T const p0[], T const p1[]) {
        T a = fabs(p0[1]);
        T b = fabs(p1[1]);
        T h = a + b;
        T t = p0[0] + a*(p1[0] - p0[0])/h;
        return t;
    }

    template<typename T>
    struct Curve {
        typedef typename T::ElementType ElementType;
        static const int LENGTH = N;
        T operator[](int i) const { return v[i]; }
        T &operator[](int i) { return v[i]; }
        T v[N];
    };
    template <typename T>
    struct vec2t {
        typedef T ElementType;
        T operator[](int i) const { return v[i]; }
        T &operator[](int i) { return v[i]; }
        T v[2];
    };
    typedef Curve<vec2t<typename ValueType::ElementType> > VECTOR;

    static int scanConvex(Real ts[], VECTOR const & p) NO_INLINE {
        typedef typename VECTOR::ElementType Scalar;
        if (!scanMinMax(p)) {
            return 0;
        }
        int n = 0;
        {
            int k = -1;
            int l = -1;
            Scalar current = std::numeric_limits<Scalar>::max();
            for (int i = 1; i < VECTOR::LENGTH; ++i) {
                if (p[i][1] * p[0][1] < 0) {
                    Scalar s = slope(&p[i][0], &p[0][0]);
                    if (s < current) {
                        current = s;
                        k = i;
                    }
                }
            }
            if (k < 0) return 0;
            current = 0;
            for (int i = 0; i < k; ++i) {
                if (p[i][1] * p[k][1] < 0) {
                    Scalar s = slope(&p[i][0], &p[k][0]);
                    if (current < s) {
                        current = s;
                        l = i;
                    }
                }
            }
            if (l < 0) return 0;
            ts[n++] = dist(&p[l][0], &p[k][0]);
        }
        {
            int k = -1;
            int l = -1;
            Scalar current = std::numeric_limits<Scalar>::max();
            for (int i = 0; i < VECTOR::LENGTH-1; ++i) {
                if (p[i][1] * p[VECTOR::LENGTH-1][1] < 0) {
                    Scalar s = slope(&p[i][0], &p[VECTOR::LENGTH-1][0]);
                    if (s < current) {
                        current = s;
                        k = i;
                    }
                }
            }
            if (k < 0) return 0;
            current = 0;
            for (int i = k+1; i < VECTOR::LENGTH; ++i) {
                if (p[i][1] * p[k][1] < 0) {
                    Scalar s = slope(&p[i][0], &p[k][0]);
                    if (current < s) {
                        current = s;
                        l = i;
                    }
                }
            }
            if (l < 0) return 0;
            ts[n++] = dist(&p[k][0], &p[l][0]);
        }
        return n;
    }

    static bool xCheck(Real rng[2], const VECTOR& curve)
    {
        Real t[2] = {0, 0};
        int nn = scanConvex(t, curve);
        if (nn) {
            if (nn != 2) return false;

            Real t0 = t[0];
            Real t1 = t[1];

            if (t0 > t1) std::swap(t0, t1);

            rng[0] = t0;
            rng[1] = t1;
            return true;
        }
        return false;
    }

    static bool getRangeU(Real out[2], PatchType const & patch) NO_INLINE {
        Real t0 = 1;
        Real t1 = 0;
        VECTOR curve;
        Real delta = 1.0/(N-1);
        for (int i = 0; i < N; ++i) {
            curve[i][0] = i*delta;
        }
        Real rng[2];
        for (int j = 0; j < N; ++j) {
            for (int i = 0; i < N; ++i) {
                curve[i][1] = (patch.Get(i,j))[1];
            }
            if (xCheck(rng, curve)) {
                if (rng[1] - rng[0] < delta) {
                    t0 = std::min(t0, rng[0]);
                    t1 = std::max(t1, rng[1]);
                } else {
                    return false;
                }
            } else {
                return false;
            }
        }
        out[0] = t0;
        out[1] = t1;
        return true;
    }

    static bool getRangeV(Real out[2], PatchType const &patch) NO_INLINE {
        Real t0 = 1;
        Real t1 = 0;
        VECTOR curve;
        Real delta = 1.0/(N-1);
        for (int i = 0; i < N; ++i) {
            curve[i][0] = i*delta;
        }
        Real rng[2];
        for (int i = 0; i < N; ++i) {
            for (int j = 0; j < N; ++j) {
                curve[j][1] = (patch.Get(i,j))[1];
            }
            if (xCheck(rng, curve)) {
                if (rng[1] - rng[0] < delta) {
                    t0 = std::min(t0,rng[0]);
                    t1 = std::max(t1,rng[1]);
                } else {
                    return false;
                }
            }else{
                return false;
            }
        }
        out[0] = t0;
        out[1] = t1;
        return true;
    }

    static bool testBezierClipU(UVT* info, PatchType const & patch,
                                Real u0, Real u1,
                                Real v0, Real v1, Real zmin, Real zmax,
                                int level, int max_level, REAL eps) NO_INLINE {
        PatchType tpatch(patch);
        rotateU(tpatch);

        ValueType min, max;
        tpatch.GetMinMax(min, max, EPS*1e-3);
        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2]) return false;//z
        //if(max[0]-min[0]<=EPS2)return false;
        //if (max[1]-min[1] <= EPS2) return false;

        //        printf("U: %f, %f\n", min[1], max[1]);

        if (isEps(min,max,eps) || isLevel(level,max_level)){
            return testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level);
        } else {
            Real tw = 1;
            Real tt0 = 1;
            Real tt1 = 0;
#if USE_BEZIERCLIP
            if (isClip(level)) {
                Real rng[2];
                if (getRangeU(rng, tpatch)) {
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }
#endif
            if (tw>=0.4) {
                PatchType tmp[2];
                patch.SplitU(tmp, 0.5);
                Real um = (u0+u1)*0.5;
                Real ut[] = {u0,um,um,u1};

                int order[2]={0,1};
                coarseSort(order, tmp);
                bool bRet = false;
                for (int i = 0; i < 2; ++i) {
                    if (testBezierClipV(info, tmp[order[i]], ut[2*order[i]], ut[2*order[i]+1], v0, v1, zmin, zmax, level+1, max_level, eps)){
                        zmax = info->t;
                        bRet = true;
                    }
                }
                return bRet;
            } else {
                tt0 = std::max<Real>(0.0, tt0-UVEPS);
                tt1 = std::min<Real>(tt1+UVEPS, 1.0);
                Real ut[] = {lerp(u0,u1,tt0),lerp(u0,u1,tt1)};
                PatchType tmp;
                patch.CropU(tmp, tt0, tt1);
                return testBezierClipV(info, tmp, ut[0], ut[1], v0, v1, zmin, zmax, level+1, max_level, eps);
            }
        }
        return false;
    }

    static bool testBezierClipV(UVT *info, PatchType const &patch,
                                Real u0, Real u1, Real v0, Real v1, Real zmin, Real zmax,
                                int level, int max_level, Real eps) NO_INLINE {
        PatchType tpatch(patch);
        rotateV(tpatch);

        ValueType min, max;
        tpatch.GetMinMax(min, max, EPS*1e-3);
        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2])return false;//z
        //if(max[0]-min[0]<=EPS2)return false;
        //if (max[1] - min[1] <= EPS2)return false;

        //        printf("V: %f, %f\n", min[1], max[1]);

        if (isEps(min,max,eps) || isLevel(level,max_level)) {
            return testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level);
        } else {
            Real tw = 1;
            Real tt0 = 1;
            Real tt1 = 0;
#if USE_BEZIERCLIP
            if (isClip(level)) {
                Real rng[2];
                if (getRangeV(rng, tpatch)) {
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }
#endif
            if (tw >= 0.4) {
                PatchType tmp[2];
                patch.SplitV(tmp, 0.5);
                Real vm = (v0+v1)*0.5;
                Real vt[] = {v0,vm,vm,v1};

                int order[2]={0,1};
                coarseSort(order, tmp);
                bool bRet = false;
                for (int i = 0; i < 2; ++i) {
                    if (testBezierClipU(info, tmp[order[i]], u0, u1,
                        vt[2*order[i]], vt[2*order[i]+1], zmin, zmax, level+1, max_level, eps)){
                        zmax = info->t;
                        bRet = true;
                    }
                }
                return bRet;
            } else {
                tt0 = std::max<Real>(0.0,tt0-UVEPS);
                tt1 = std::min<Real>(tt1+UVEPS,1.0);
                Real vt[] = {lerp(v0,v1,tt0),lerp(v0,v1,tt1)};
                PatchType tmp;
                patch.CropV(tmp, tt0, tt1);
                return testBezierClipU(info, tmp, u0, u1, vt[0], vt[1], zmin, zmax, level+1, max_level, eps);
            }
        }
        return false;
    }

    static bool testBezierClipL(UVT* info, PatchType const &patch,
                                Real u0, Real u1, Real v0, Real v1, Real zmin, Real zmax, int level) NO_INLINE {
        Real t = Real(0), u = Real(0), v = Real(0);
#if DIRECT_BILINEAR
        ValueType P[4];
        P[0] = patch.Get(0, 0);
        P[1] = patch.Get(N-1, 0);
        P[2] = patch.Get(0, N-1);
        P[3] = patch.Get(N-1, N-1);
        if (testBilinearPatch(&t, &u, &v, P, zmin, zmax)) {
            u = u0*(1-u)+u1*u;
            v = v0*(1-v)+v1*v;
            info->u = u;
            info->v = v;
            info->t = t;
            info->level = level;
            return true;
        }
        return false;
#else
        bool bRet = false;
        int nPu = N - 1;
        int nPv = N - 1;

        REAL du = REAL(1)/nPu;
        REAL dv = REAL(1)/nPv;

        REAL uu = 0;
        REAL vv = 0;

        ValueType P[4];
        for (int j = 0; j < nPv; ++j) {
            for (int i = 0; i < nPu; ++i) {
                P[0] = patch.Get(i  ,j  );
                P[1] = patch.Get(i+1,j  );
                P[2] = patch.Get(i  ,j+1);
                P[3] = patch.Get(i+1,j+1);
                if(testBilinearPatch(&t, &u, &v, P, zmin, zmax)){

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
                    info->level = level;
                    bRet = true;
                }
            }
        }

        if(bRet) {
            ValueType p = patch.Evaluate(uu,vv);
            info->t = p[2];
        }

        return bRet;
#endif
    }
    bool intersectAABB(RangeAABB *rng, ValueType const & min, ValueType const & max,
                       Ray const & r, Real tmin, Real tmax) const NO_INLINE {

        //int phase = r.phase();
        int sign[3] = {r.dirSign[0],r.dirSign[1],r.dirSign[2]};
        ValueType box[2] = { min, max };

        ValueType org(r.org[0], r.org[1], r.org[2]);
        ValueType idir(r.invDir[0], r.invDir[1], r.invDir[2]);

        for (int i = 0; i < 3; ++i) {
            tmin = std::max(tmin, (box[  sign[i]][i]-org[i])*idir[i]);
            tmax = std::min(tmax, (box[1-sign[i]][i]-org[i])*idir[i]);
        }
        rng->tmin = tmin;
        rng->tmax = tmax;

        return rng->tmin <= rng->tmax;
    }


    PatchType _patch;
    Real _uRange[2];
    Real _vRange[2];
    ValueType _min;
    ValueType _max;
    Real _eps;
    int _maxLevel;
};

}   // end OsdUtil

//}  // end namespace OPENSUBDIV_VERSION
//using namespace OPENSUBDIV_VERSION;

//}  // end namespace OpenSubdiv

#endif  // OSDUTIL_BEZIER_INTERSECT_H
