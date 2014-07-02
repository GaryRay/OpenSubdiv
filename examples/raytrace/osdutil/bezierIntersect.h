#ifndef OSDUTIL_BEZIER_INTERSECT_H
#define OSDUTIL_BEZIER_INTERSECT_H

#include "bezier.h"
#include "bilinearIntersect.h"

#include <limits>
#include <cfloat>
#include "../common.h"

//namespace OpenSubdiv {
//namespace OPENSUBDIV_VERSION {

#define USE_BEZIERCLIP 1
#define USE_COARSESORT 1
#define USE_UVCROP 1

namespace OsdUtil {

template <typename REAL>
struct Epsilon {
    static const REAL EPS = 1e-4;
    static const REAL UVEPS = 1.0/32.0;
};
template <>
struct Epsilon<double> {
    static const double EPS = 1e-4;    // Emprically it seems work(@syoyo). At least 1e-6 or smaller cause some crack. 2014/06/15. revisit
    static const double UVEPS = 1.0/32.0;
};

#define Get16Bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
                       +(uint32_t)(((const uint8_t *)(d))[0]) )

inline uint32_t fastHash(const char *data, int len, uint32_t hash) {
    uint32_t tmp;
    int rem;

    rem = len & 3;
    len >>= 2;

    /* Main loop */
    for (;len > 0; len--) {
        hash  += Get16Bits (data);
        tmp    = (Get16Bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }

    /* Handle end cases */
    switch (rem) {
    case 3: hash += Get16Bits (data);
        hash ^= hash << 16;
        hash ^= data[sizeof (uint16_t)] << 18;
        hash += hash >> 11;
        break;
    case 2: hash += Get16Bits (data);
        hash ^= hash << 11;
        hash += hash >> 17;
        break;
    case 1: hash += *data;
        hash ^= hash << 10;
        hash += hash >> 1;
    }

    /* Force "avalanching" of final 127 bits */
    hash ^= hash << 3;
    hash += hash >> 5;
    hash ^= hash << 4;
    hash += hash >> 17;
    hash ^= hash << 25;
    hash += hash >> 6;

    return hash;
}

template<class VALUE_TYPE, class REAL, int N=4>
class OsdUtilBezierPatchIntersection {
public:

    typedef VALUE_TYPE ValueType;
    typedef REAL Real;
    typedef OsdUtilBezierPatch<ValueType, Real, N> PatchType;

    static const REAL EPS = Epsilon<Real>::EPS;
    static const REAL UVEPS = Epsilon<Real>::UVEPS;
    static const int  DEFAULT_MAX_LEVEL = 10;

    struct RangeAABB {
        Real tmin, tmax;
    };
    struct UVT {
        Real u, v, t;
        int level;
        int quadHash;
    };

    OsdUtilBezierPatchIntersection(PatchType const &patch) NO_INLINE :
        _patch(patch),
        _eps(EPS),
        _maxLevel(DEFAULT_MAX_LEVEL),
        _uvMargin(true),
        _cropUV(true), _useBezierClip(true),
        _useTriangle(false),
        _directBilinear(false),
        _wcpFlag(0)
    {

        _uRange[0] = _vRange[0] = 0;
        _uRange[1] = _vRange[1] = 1;

        patch.GetMinMax(_min, _max, Real(1e-3));
    }
    ~OsdUtilBezierPatchIntersection() { }

    void SetEpsilon(REAL eps){
        _eps = std::max(std::numeric_limits<REAL>::epsilon()*REAL(1e+1), eps);
    }
    void SetMaxLevel(int maxLevel){
        _maxLevel = maxLevel;
    }
    void SetUVMergin(REAL uvMargin){
        _uvMargin = uvMargin;
    }
    void SetCropUV(bool cropUV){
        _cropUV = cropUV;
    }
    void SetUseBezierClip(bool useBezierClip){
        _useBezierClip = useBezierClip;
    }
    void SetUseTriangle(bool useTriangle){
        _useTriangle = useTriangle;
    }
    void SetDirectBilinear(bool directBilinear){
        _directBilinear = directBilinear;
    }
    void SetWatertightFlag(int wcpFlag) {
        _wcpFlag = wcpFlag;
    }

    bool Test(Intersection* info, Ray const &r, Real tmin, Real tmax) NO_INLINE {

        RangeAABB rng;
        if (intersectAABB(&rng, _min, _max, r, tmin, tmax)) {
            tmin = std::max(tmin, rng.tmin);
            tmax = std::min(tmax, rng.tmax);

            typename ValueType::Matrix4Type mat(ValueType(r.org), ValueType(r.dir)); // getZAlign
            _mat = mat;

            return testInternal(info, r, tmin, tmax);
        }
        return false;
    }

    Real ComputeEpsilon(Ray const & r, Real eps)const
    {
        RangeAABB rng;
        if(r.hasDifferential && intersectAABB(&rng, _min, _max, r, 0, std::numeric_limits<Real>::max()))
        {
            Real t = std::max(Real(0), rng.tmin);
            //ValueType org(r.org[0], r.org[1], r.org[2]);
            ValueType dir(r.dir[0], r.dir[1], r.dir[2]);
            ValueType dirDX = (dir + ValueType(r.dDdx[0],r.dDdx[1],r.dDdx[2]));
            ValueType dirDY = (dir + ValueType(r.dDdy[0],r.dDdy[1],r.dDdy[2]));
            dirDX.normalize();
            dirDY.normalize();
            ValueType P0 = t * dir;
            ValueType PX = t * dirDX;
            ValueType PY = t * dirDY;

            Real lx = (PX-P0).length();
            Real ly = (PY-P0).length();

            eps = std::min(lx,ly)*Real(0.25);
        }
        return eps;
    }

protected:
    bool testInternal(Intersection* info, const Ray& r, Real tmin, Real tmax) const NO_INLINE {
        UVT uvt;
        PatchType patch(_patch);
        patch.Transform(_mat);
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
            info->quadHash = uvt.quadHash;
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
            trace("hit t = %f, uv = (%f, %f)\n", t, u, v);

            return true;
        }
        return false;
    }
    bool testBezierPatch(UVT* info, PatchType const & patch, Real zmin, Real zmax, Real eps, bool wcpPass=false) const {
        ValueType min, max;
        patch.GetMinMax(min, max, eps*1e-3);

        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2]) return false;//z

        if (_cropUV) return testBezierClipRangeU(info, patch, 0, 1, 0, 1, zmin, zmax, 0, _maxLevel, eps, wcpPass);
        else         return testBezierClipU     (info, patch, 0, 1, 0, 1, zmin, zmax, 0, _maxLevel, eps, wcpPass);
    }

    // ----------------------------------------------------------------------
    static void rotateU(PatchType& patch) NO_INLINE
    {
        ValueType dx = patch.GetLv();
        typename ValueType::Matrix3Type rot(dx);
        patch.Transform(rot);
    }

    static void rotateV(PatchType& patch) NO_INLINE
    {
        ValueType dx = patch.GetLu();
        typename ValueType::Matrix3Type rot(dx);
        patch.Transform(rot);
    }

    static bool isEps(ValueType const & min, ValueType const & max, Real eps) {
        //return ((max[1]-min[1])<=eps);
        REAL xw = max[0]-min[0];
        REAL yw = max[1]-min[1];
        return (xw<=eps && yw<=eps);
    }

    static bool isLevel(int level, int max_level) {
        return (level>=max_level);
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
    static Real lerp(Real a, Real b, Real u) {
        Real t = Real(1) - (Real(1)-u);
        Real s = Real(1) - t;
        return s*a+t*b;
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
        T d1 = b[1] - a[1];
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
            if (k >= 0){
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
                if (l >= 0){
                    ts[n++] = dist(&p[l][0], &p[k][0]);
                }
            }
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
            if (k >= 0){
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
                if (l >= 0) {
                    ts[n++] = dist(&p[k][0], &p[l][0]);
                }
            }
        }
        return n;
    }

    static bool xCheck(Real rng[2], const VECTOR& curve)
    {
        Real t[4] = {0};
        int nn = scanConvex(t, curve);
        if (nn) {
            Real t0 = t[0];
            Real t1 = t[0];

            for(int i=1;i<nn;i++){
                t0 = std::min(t0, t[i]);
                t1 = std::max(t1, t[i]);
            }

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
                t0 = std::min(t0, rng[0]);
                t1 = std::max(t1, rng[1]);
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
                t0 = std::min(t0,rng[0]);
                t1 = std::max(t1,rng[1]);
            }else{
                return false;
            }
        }
        out[0] = t0;
        out[1] = t1;
        return true;
    }

    bool testBezierClipU(UVT* info, PatchType const & patch,
                         Real u0, Real u1, Real v0, Real v1,
                         Real zmin, Real zmax,
                         int level, int max_level, Real eps, bool wcpPass) const NO_INLINE {
        PatchType tpatch(patch);
        rotateU(tpatch);

        trace("testBezierClipU (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
              u0, u1, v0, v1, zmin, zmax, level);

        ValueType min, max;
        tpatch.GetMinMax(min, max, eps*1e-3);
        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2]) return false;//z
        
        bool bClip = isClip(level);
        if (bClip && (isEps(min,max,eps) || isLevel(level,max_level))){
            return testBezierClipL2(info, patch, u0, u1, v0, v1, zmin, zmax, level, wcpPass);
        } else {
            Real tw = 1;
            Real tt0 = 0;
            Real tt1 = 1;

            if (_useBezierClip & bClip) {
                Real rng[2];
                if (getRangeU(rng, tpatch)) {
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }

            if (tw>=0.4) {
                PatchType tmp[2];
                patch.SplitU(tmp, 0.5);
                Real um = (u0+u1)*0.5;
                Real ut[] = {u0,um,um,u1};

                int order[2]={0,1};
                coarseSort(order, tmp);
                bool bRet = false;
                for (int i = 0; i < 2; ++i) {
                    if (testBezierClipV(info, tmp[order[i]], ut[2*order[i]], ut[2*order[i]+1],
                                        v0, v1, zmin, zmax, level+1, max_level, eps, wcpPass)){
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
                return testBezierClipV(info, tmp, ut[0], ut[1],
                                       v0, v1, zmin, zmax, level+1, max_level, eps, wcpPass);
            }
        }
        return false;
    }

    bool testBezierClipV(UVT *info, PatchType const &patch,
                         Real u0, Real u1, Real v0, Real v1, Real zmin, Real zmax,
                         int level, int max_level, Real eps, bool wcpPass) const NO_INLINE {
        PatchType tpatch(patch);
        rotateV(tpatch);

        trace("testBezierClipV (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
              u0, u1, v0, v1, zmin, zmax, level);

        ValueType min, max;
        tpatch.GetMinMax(min, max, eps*1e-3);
        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2])return false;//z

        bool bClip = isClip(level);
        if (bClip && (isEps(min,max,eps) || isLevel(level,max_level))) {
            return testBezierClipL2(info, patch, u0, u1, v0, v1, zmin, zmax, level, wcpPass);
        } else {
            Real tw = 1;
            Real tt0 = 0;
            Real tt1 = 1;

            if (_useBezierClip & bClip) {
                Real rng[2];
                if (getRangeV(rng, tpatch)) {
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }

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
                                        vt[2*order[i]], vt[2*order[i]+1],
                                        zmin, zmax, level+1, max_level, eps, wcpPass)){
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
                return testBezierClipU(info, tmp, u0, u1,
                                       vt[0], vt[1], zmin, zmax, level+1, max_level, eps, wcpPass);
            }
        }
        return false;
    }

    bool testBezierClipL2(UVT* info, PatchType const &patch,
                          Real u0, Real u1, Real v0, Real v1,
                          Real zmin, Real zmax, int level, bool wcpPass) const NO_INLINE {
        
        if (testBezierClipL(info, patch, u0, u1, v0, v1, zmin, zmax, level)) {
            return true;
        }
        
        if (_wcpFlag != 0 && !wcpPass && (u0 == 0 || u1 ==1 || v0 == 0 || v1 == 1)) {
            // TODO: more efficient test (using wcpFlag and uv)

            // test against split faces too
                /*
                  +---+---+
               ^  | 2 | 3 |
               |  +---+---+
               |  | 0 | 1 |
               u  +---+---+
                    v---->
                 */
            PatchType tmp[2];
            PatchType children[4];
            _patch.SplitU(tmp, 0.5);
            tmp[0].SplitV(&children[0], 0.5);
            tmp[1].SplitV(&children[2], 0.5);
            children[0].Transform(_mat);
            children[1].Transform(_mat);
            children[2].Transform(_mat);
            children[3].Transform(_mat);
            Real uvs[4][4] = { { u0*2, u1*2, v0*2, v1*2},
                               { u0*2, u1*2, (v0-0.5)*2, (v1-0.5)*2 },
                               { (u0-0.5)*2, (u1-0.5)*2, v0*2, v1*2 },
                               { (u0-0.5)*2, (u1-0.5)*2, (v0-0.5)*2, (v1-0.5)*2 } };
            PatchType cp;
            trace("Subface uv (%f-%f),(%f-%f)\n", u0, u1, v0, v1);

            // XXX: revisit this logic... seems very redundant.
            
            for (int i = 0; i < 4; ++i) {
                if (uvs[i][0] > 1 || uvs[i][1] < 0 ||
                    uvs[i][2] > 1 || uvs[i][3] < 0) continue;
                if (u0 == 0 && (i == 2 || i == 3)) continue;
                if (u1 == 1 && (i == 0 || i == 1)) continue;
                if (v0 == 0 && (i == 1 || i == 3)) continue;
                if (v1 == 1 && (i == 0 || i == 2)) continue;

                trace("Child pass %d\n", i);
                if (testBezierPatch(info, children[i], zmin, zmax, _eps, /*wcpPass=*/true)) {
                    if (i == 0) {
                        info->u = info->u * 0.5;
                        info->v = info->v * 0.5;
                    } else if (i == 1) {
                        info->u = info->u * 0.5;
                        info->v = info->v * 0.5 + 0.5;
                    } else if (i == 2) {
                        info->u = info->u * 0.5 + 0.5;
                        info->v = info->v * 0.5;
                    } else {
                        info->u = info->u * 0.5 + 0.5;
                        info->v = info->v * 0.5 + 0.5;
                    }
                    return true;
                }
            }
        }
        return false;
    }

    bool testBezierClipL(UVT* info, PatchType const &patch,
                         Real u0, Real u1, Real v0, Real v1,
                         Real zmin, Real zmax, int level) const NO_INLINE {
        Real t = Real(0), u = Real(0), v = Real(0);

        trace("testBezierClipL (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
              u0, u1, v0, v1, zmin, zmax, level);

        if (_directBilinear) {
            ValueType P[4];
            P[0] = patch.Get(0, 0);
            P[1] = patch.Get(N-1, 0);
            P[2] = patch.Get(0, N-1);
            P[3] = patch.Get(N-1, N-1);
            if (!_useTriangle) {
                if (testBilinearPatch(&t, &u, &v, P, zmin, zmax, _uvMargin)) {
                    u = u0*(1-u)+u1*u;
                    v = v0*(1-v)+v1*v;
                    info->u = u;
                    info->v = v;
                    info->t = t;
                    info->level = level;
                    info->quadHash = computeHash(u0, u1, v0, v1);
                    return true;
                }
            } else {
                if (testQuadPlane(&t, &u, &v, P, ValueType(0,0,0), ValueType(0,0,1), zmin, zmax)){
                    u = u0*(1-u)+u1*u;
                    v = v0*(1-v)+v1*v;
                    info->u = u;
                    info->v = v;
                    info->t = t;
                    info->level = level;
                    info->quadHash = computeHash(u0, u1, v0, v1);
                    return true;
                }
            }
            return false;
        }

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
                if(!_useTriangle){
                    if (testBilinearPatch(&t, &u, &v, P, zmin, zmax, _uvMargin)){

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
                        info->quadHash = computeHash(u0, u1, v0, v1);
                        bRet = true;
                    }
                }else{
                    if (testQuadPlane(&t, &u, &v, P, ValueType(0,0,0), ValueType(0,0,1), zmin, zmax)){
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
                        info->quadHash = computeHash(u0, u1, v0, v1);
                        bRet = true;
                    }
                }
            }
        }

        if(bRet) {
            ValueType p = patch.Evaluate(uu,vv);
            info->t = p[2];
        }

        return bRet;
    }

    bool testBezierClipRangeU(UVT* info, PatchType const & patch,
                              Real u0, Real u1,
                              Real v0, Real v1, Real zmin, Real zmax,
                              int level, int max_level, Real eps, bool wcpPass) const NO_INLINE {
        PatchType mpatch(patch, u0, u1, v0, v1);
        PatchType tpatch(mpatch);
        rotateU(tpatch);

        trace("testBezierClipRangeU (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
              u0, u1, v0, v1, zmin, zmax, level);

        ValueType min, max;
        tpatch.GetMinMax(min, max, eps*1e-3);
        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2]) return false;//z
        
        bool bClip = isClip(level);
        if (bClip && (isEps(min,max,eps) || isLevel(level,max_level))){
            return testBezierClipL2(info, mpatch, u0, u1, v0, v1, zmin, zmax, level, wcpPass);
        } else {
            Real tw = 1;
            Real tt0 = 0;
            Real tt1 = 1;

            if (_useBezierClip & bClip) {
                Real rng[2];
                if (getRangeU(rng, tpatch)) {
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }

            if (tw>=0.4) {
                Real um = (u0+u1)*0.5;
                um = Real(1.0) - (Real(1.0) - um);
                Real ut[] = {u0,um,um,u1};
                int order[2]={0,1};
                bool bRet = false;
                for (int i = 0; i < 2; ++i) {
                    if (testBezierClipRangeV(info, patch, ut[2*order[i]], ut[2*order[i]+1],
                                             v0, v1, zmin, zmax, level+1, max_level, eps, wcpPass)) {
                        zmax = info->t;
                        bRet = true;
                    }
                }
                return bRet;
            } else {
                tt0 = std::max<Real>(0.0, tt0-UVEPS);
                tt1 = std::min<Real>(tt1+UVEPS, 1.0);
                Real ut[] = {lerp(u0,u1,tt0),lerp(u0,u1,tt1)};
                return testBezierClipRangeV(info, patch, ut[0], ut[1],
                                            v0, v1, zmin, zmax, level+1, max_level, eps, wcpPass);
            }
        }
        return false;
    }

    bool testBezierClipRangeV(UVT *info, PatchType const &patch,
                              Real u0, Real u1, Real v0, Real v1, Real zmin, Real zmax,
                              int level, int max_level, Real eps, bool wcpPass) const NO_INLINE {
        PatchType mpatch(patch, u0, u1, v0, v1);
        PatchType tpatch(mpatch);
        rotateV(tpatch);

        trace("testBezierClipRangeV (%f, %f) - (%f, %f) z:%f, %f  level=%d\n",
              u0, u1, v0, v1, zmin, zmax, level);

        ValueType min, max;
        tpatch.GetMinMax(min, max, eps*1e-3);
        if (0 < min[0] || max[0] < 0) return false;//x
        if (0 < min[1] || max[1] < 0) return false;//y
        if (max[2] < zmin || zmax < min[2])return false;//z

        bool bClip = isClip(level);
        if (bClip && (isEps(min,max,eps) || isLevel(level,max_level))) {
            return testBezierClipL2(info, mpatch, u0, u1, v0, v1, zmin, zmax, level, wcpPass);
        } else {
            Real tw = 1;
            Real tt0 = 0;
            Real tt1 = 1;

            if (_useBezierClip & bClip) {
                Real rng[2];
                if (getRangeV(rng, tpatch)) {
                    tt0 = rng[0];
                    tt1 = rng[1];
                    tw = tt1-tt0;
                }
            }

            if (tw >= 0.4) {
                Real vm = (v0+v1)*0.5;
                vm = Real(1.0) - (Real(1.0) - vm);
                Real vt[] = {v0,vm,vm,v1};
                int order[2]={0,1};
                bool bRet = false;
                for (int i = 0; i < 2; ++i) {
                    if (testBezierClipRangeU(info, patch, u0, u1,
                                             vt[2*order[i]], vt[2*order[i]+1],
                                             zmin, zmax, level+1, max_level, eps, wcpPass)){
                        zmax = info->t;
                        bRet = true;
                    }
                }
                return bRet;
            } else {
                tt0 = std::max<Real>(0.0,tt0-UVEPS);
                tt1 = std::min<Real>(tt1+UVEPS,1.0);
                Real vt[] = {lerp(v0,v1,tt0),lerp(v0,v1,tt1)};
                return testBezierClipRangeU(info, patch, u0, u1,
                                            vt[0], vt[1], zmin, zmax, level+1, max_level, eps, wcpPass);
            }
        }
        return false;
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

    template<typename T>
    uint32_t computeHash(T a, T b, T c, T d) const {
        uint32_t hash = fastHash((const char *)&_patch, sizeof(_patch), 0);
        T v[4] = {a, b, c, d};
        hash = fastHash((const char*)v, sizeof(v), hash);
        return hash;
    }

    PatchType _patch;
    Real _uRange[2];
    Real _vRange[2];
    ValueType _min;
    ValueType _max;
    Real _eps;
    int _maxLevel;
    float _uvMargin;
    bool _cropUV;
    bool _useBezierClip;
    bool _useTriangle;
    bool _directBilinear;
    int _wcpFlag;
    typename ValueType::Matrix4Type _mat;
};

}   // end OsdUtil

//}  // end namespace OPENSUBDIV_VERSION
//using namespace OPENSUBDIV_VERSION;

//}  // end namespace OpenSubdiv

#endif  // OSDUTIL_BEZIER_INTERSECT_H
