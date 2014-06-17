//
// Ray differential utils.
// Grabbed from lucille and pbrt.
//

#include "raydiff.h"

namespace {

const real kEPS = 1.0e-4;

bool
SolveLinearSystem22(const real A[2][2], const real B[2], real x[2])
{
    real det = A[0][0]*A[1][1] - A[0][1]*A[1][0];
    if (std::abs(det) < kEPS)
        return false;
    real invDet = 1.0/det;
    x[0] = (A[1][1]*B[0] - A[0][1]*B[1]) * invDet;
    x[1] = (A[0][0]*B[1] - A[1][0]*B[0]) * invDet;
    return true;

}

//
// For a input triangle vertex p0, p1 and p2, compute plane L0, L1 and L2
// whose properties are:
//
// L0 : p2 . L0 = 1, L0 contains p0 and p1, and perp. to the triangle.
// L1 : p0 . L1 = 1, L1 contains p1 and p2, and perp. to the triangle.
// L2 : p1 . L2 = 1, L2 contains p2 and p0, and perp. to the triangle.
//
void CalculatePlaneCoefficients(
  real3& L0,
  real3& L1,
  real3& L2,
  const real3 v0,
  const real3 v1,
  const real3 v2)
{
  real3 e1 = v1 - v0;
  real3 e2 = v2 - v0;

  real3 e10 = v1 - v0;
  real3 e21 = v2 - v1;
  real3 e02 = v0 - v2;

  // Triangle normal
  real3 n = vcross(e1, e2);
  n.normalize();

  // Normal of L0 = n x e21
  // Normal of L1 = n x e02
  // Normal of L2 = n x e10
  real3 Ln0 = vcross(n, e21); Ln0.normalize();
  real3 Ln1 = vcross(n, e02); Ln1.normalize();
  real3 Ln2 = vcross(n, e10); Ln2.normalize();

    // To make (L0 . p0) + d = 1, calculate normalization factor w0 for
    // plane coefficients of L0 which is computed as
    // w0 = 1.0 / ((L0 . p0) + d)
    //
    // Normalization coeff for L1 and L2 are computed in similar fashion.

    real d0 = -(vdot(v2, Ln0));
    real d1 = -(vdot(v0, Ln1));
    real d2 = -(vdot(v1, Ln2));

    real w0_denom = vdot(Ln0, v0) + d0;
    real w1_denom = vdot(Ln1, v1) + d1;
    real w2_denom = vdot(Ln2, v2) + d2;

    real w0, w1, w2;

    if (fabs(w0_denom) > kEPS) {
        w0 = 1.0 / w0_denom;
    } else {
        w0 = 1.0;
    }

    if (fabs(w1_denom) > kEPS) {
        w1 = 1.0 / w1_denom;
    } else {
        w1 = 1.0;
    }

    if (fabs(w2_denom) > kEPS) {
        w2 = 1.0 / w2_denom;
    } else {
        w2 = 1.0;
    }

    L0 = Ln0 * w0;
    L1 = Ln1 * w1;
    L2 = Ln2 * w2;

    L0[3] = d0 * w0;
    L1[3] = d1 * w1;
    L2[3] = d2 * w2;




#if 0
   // Validate
   {
        real dl0 = L0.dot(v2) + L0[3];
        real dl1 = L1.dot(v0) + L1[3];
        real dl2 = L2.dot(v1) + L2[3];

        //assert(fabs(dl0) < (10.0 * kEPS));
        //assert(fabs(dl1) < (10.0 * kEPS));
        //assert(fabs(dl2) < (10.0 * kEPS));

        real dot0 = L0.dot(v0) + d0 * w0;
        real dot1 = L1.dot(v1) + d1 * w1;
        real dot2 = L2.dot(v2) + d2 * w2;

        //assert(fabs(fabs(dot0) - 1.0) < (10.0 * kEPS));
        //assert(fabs(fabs(dot1) - 1.0) < (10.0 * kEPS));
        //assert(fabs(fabs(dot2) - 1.0) < (10.0 * kEPS));
    }
#endif

}

void
CalculateDifferencesForInterpolatedNormal(
    real3       &dNdx,      /// [out]
    real3       &dNdy,      /// [out]
    const real3 &P,
    const real3 &dPdx,
    const real3 &dPdy,
    const real3 &n0,
    const real3 &n1,
    const real3 &n2,
    const real3 &L0,
    const real3 &L1,
    const real3 &L2)
{
    // n = (L0 . P) N0 + (L1 . P) N1 + (N2 . P) N2
    //
    //           n
    // N = -------------
    //     (n . n)^(1/2)
    //
    //
    // dn         dP             dP             dP
    // -- = (L0 . --) N0 + (L1 . --) N1 + (L2 . --) N2
    // dx         dx             dx             dx
    //
    // dN   (n.n) (dn/dx) - (n . (dn/dx)) n
    // -- = -------------------------------
    // dx             (n.n)^(3/2)

    real dotL0P = vdot(L0, P) + L0[3];    // homogeneous
    real dotL1P = vdot(L1, P) + L1[3];
    real dotL2P = vdot(L2, P) + L2[3];

    real3 n = dotL0P * n0 + dotL1P * n1 + dotL2P * n2;

    // printf("calc(n) = %f, %f, %f\n", n[0], n[1], n[2]);

    real dotL0dPdx = vdot(L0, dPdx);
    real dotL1dPdx = vdot(L1, dPdx);
    real dotL2dPdx = vdot(L2, dPdx);

    real dotL0dPdy = vdot(L0, dPdy);
    real dotL1dPdy = vdot(L1, dPdy);
    real dotL2dPdy = vdot(L2, dPdy);

    //printf("L0.P = %f\n", dotL0dPdx);
    //printf("L1.P = %f\n", dotL1dPdx);
    //printf("L2.P = %f\n", dotL2dPdx);

    real3 dndx = dotL0dPdx * n0 + dotL1dPdx * n1 + dotL2dPdx * n2;
    real3 dndy = dotL0dPdy * n0 + dotL1dPdy * n1 + dotL2dPdy * n2;

    real dotnn = vdot(n, n);
    real dotndndx = vdot(n, dndx);
    real dotndndy = vdot(n, dndy);

    // dN   (n.n) (dn/dx) - (n . (dn/dx)) n
    // -- = -------------------------------
    // dx             (n.n)^(3/2)

    dNdx = (dotnn * dndx - dotndndx * n) / pow(dotnn, (3.0 / 2.0));
    dNdy = (dotnn * dndy - dotndndy * n) / pow(dotnn, (3.0 / 2.0));

    //printf("dNdx = %f, %f, %f\n", dNdx[0], dNdx[1], dNdx[2]);
    //printf("dNdy = %f, %f, %f\n", dNdy[0], dNdy[1], dNdy[2]);
    
}

inline void
FindBarycentricCoord(
    real &u,
    real &v,
    real &t,
    const real3 &v0,
    const real3 &v1,
    const real3 &v2,
    const real3 &rayorg,
    const real3 &raydir)
{
    real3  e1, e2;
    real3  p, s, q;
    real       a, inva = 0.0;

    e1 = v1 - v0;
    e2 = v2 - v0;

    p = vcross(raydir, e2);

    a = vdot(e1, p);

    if (std::abs(a) > kEPS) {
        inva = 1.0 / a;
    }

    (s = rayorg - v0);
    q = vcross(s, e1);

    u = vdot(s, p) * inva;
    v = vdot(q, raydir) * inva;
    t = vdot(e2, q) * inva;
}

} // namespace

using namespace mallie;

bool
ComputeDifferentialsForReflection(
    real3       &dDdx,                
    real3       &dDdy,                
    const Ray   &ray,
    const real3 &P,
    const real3 &D,
    const real3 &N,
    const real3 &p0,
    const real3 &p1,
    const real3 &p2,
    const real3 &n0,
    const real3 &n1,
    const real3 &n2)
{
    //
    // dP'   dP
    // --  = --
    // dx    dx
    //
    // dD'   dD     |       dN   d(D.N)   | 
    // --  = -- - 2 | (D.N) -- + ------ N |
    // dx    dx     |       dx     dx     |
    //
    // where
    //
    //   d(D.N)   dD        dN
    //   ------ = -- N  + D --
    //     dx     dx        dx

    if (!ray.hasDifferential) {

        dDdx[0] = dDdx[1] = dDdx[2] = 0.0;
        dDdy[0] = dDdy[1] = dDdy[2] = 0.0;

        return false;

    } else {

        real3   L0, L1, L2;

        CalculatePlaneCoefficients(
            L0, L1, L2,
            p0, p1, p2);

        real3 dNdx, dNdy;

        real3 dPdx = ray.rx.dP;
        real3 dPdy = ray.ry.dP;
        real3 dD_dx = ray.rx.dD;
        real3 dD_dy = ray.ry.dD;

        CalculateDifferencesForInterpolatedNormal(
            dNdx, dNdy,
            P,
            dPdx, dPdy,
            n0, n1, n2, L0, L1, L2);
            
        real3 dDNdx = dDdx * N + D * dNdx;
        real3 dDNdy = dDdy * N + D * dNdy;

        real dotDN = vdot(D, N);

        //real3 R = D - 2.0 * dotDN * N;
        // printf("reflect vec = %f, %f, %f\n", R[0], R[1], R[2]);

        dDdx = dD_dx - (2.0 * (dotDN * dNdx + dDNdx * N));
        dDdy = dD_dy - (2.0 * (dotDN * dNdy + dDNdy * N));

    }

    return true;

}

bool
ComputeDifferentialsForRefraction(
    real3       &dDdx,                
    real3       &dDdy,                
    const Ray       &ray,
    const real3 &P,
    const real3 &D,
    const real3 &N,
    const real3 &p0,
    const real3 &p1,
    const real3 &p2,
    const real3 &n0,
    const real3 &n1,
    const real3 &n2,
    real             eta)
{
    // P' = P
    // D' = eta D - mu N 
    //
    // mu = eta ((D . N) - (D' . N))
    //
    // D' . N = -sqrt(1 - eta^2 (1 - (D.N)^2))
    //
    //
    // dP'   dP
    // --  = --
    // dx    dx
    //
    // dD'       dD   |    dN   dmu   | 
    // --  = eta -- - | mu -- + --- N |
    // dx        dx   |    dx   dx    |
    //
    // where
    //
    //   dmu   |       eta^2 (D.N) | d(D.N)
    //   --- = | eta - ----------- | ------
    //   dx    |         (D'.N)    |   dx

    if (!ray.hasDifferential) {

        dDdx[0] = dDdx[1] = dDdx[2] = 0.0;
        dDdy[0] = dDdy[1] = dDdy[2] = 0.0;

        return false;

    } else {

        real3   L0, L1, L2;

        CalculatePlaneCoefficients(
            L0, L1, L2,
            p0, p1, p2);

        real3 dNdx, dNdy;
        real3 dPdx, dPdy;

        dPdx = ray.rx.dP;
        dPdy = ray.ry.dP;

        // dN
        // --
        // dx
        CalculateDifferencesForInterpolatedNormal(
            dNdx, dNdy,
            P,
            dPdx, dPdy,
            n0, n1, n2, L0, L1, L2);
            
        real dotDN = vdot(D, N);
        real3 NN = N;

        if (dotDN < 0.0) {
            dotDN = -dotDN;
            eta = 1.0 / eta;
        } else {
            NN = N.neg();
        }

        //   d(D.N)   dD        dN
        //   ------ = -- N  + D --
        //     dx     dx        dx
        real3 dDNdx = dDdx * NN + dNdx * D;
        real3 dDNdy = dDdy * NN + dNdy * D;

        // dotRN = D'.N
        real k = 1.0 - eta * eta * (1.0 - dotDN * dotDN);
        real dotRN = 0.0;

        // printf("k = %f\n", k);

        if (k > 0.0) {
            dotRN = sqrt(k);
        } else {
            // @todo { Total Internal Reflection, calculate diff for reflection ray. }
        }

        real mu = eta * dotDN - dotRN;

        //real3 R = eta * D + mu * N;  // D'

        //printf("Diff: eta = %f\n", eta);
        //printf("Diff: D.N = %f\n", dotDN);
        //printf("Diff: D = %f, %f, %f\n", D[0], D[1], D[2]);
        //printf("Diff: N = %f, %f, %f\n", N[0], N[1], N[2]);
        //printf("Diff: refract vec = %f, %f, %f\n", R[0], R[1], R[2]);

        //   dmu   |       eta^2 (D.N) | d(D.N)
        //   --- = | eta - ----------- | ------
        //   dx    |         (D'.N)    |   dx

        real3 dmudx, dmudy;

        dmudx = (eta - (eta*eta*dotDN) / dotRN) * dDNdx;
        dmudy = (eta - (eta*eta*dotDN) / dotRN) * dDNdy;

        // dD'       dD   |    dN   dmu   | 
        // --  = eta -- - | mu -- + --- N |
        // dx        dx   |    dx   dx    |

        real3 dD_dx = ray.rx.dD;
        real3 dD_dy = ray.ry.dD;

        dDdx = eta * dD_dx - (mu * dNdx + dmudx * N); 
        dDdy = eta * dD_dy - (mu * dNdy + dmudy * N); 

    }

    return true;

}

bool
ComputeDifferentialsForTransfer(
    real3 &dPdx,                
    real3 &dPdy,                
    real      &dudx,                 
    real      &dvdx,                 
    real      &dudy,                 
    real      &dvdy,                 
    const Ray &ray,
    const real3 &p,
    const real3 &n,
    const real3 &dPdu,          ///< dP/du (tangent)
    const real3 &dPdv)          ///< dP/dv (binormal)
{
    if (!ray.hasDifferential) {

        dPdx[0] = dPdx[1] = dPdx[2] = 0.0;
        dPdy[0] = dPdy[1] = dPdy[2] = 0.0;
        dudx = dvdx = 0.0;
        dudy = dvdy = 0.0;

        return false;

    } else {

        // Compute aux points px and py.
        real d = -(vdot(n, p));
        real3 rxv = ray.rx.P;
        real tx = -(vdot(n, rxv) + d) / vdot(n, ray.rx.D);
        real3 px = ray.rx.P + tx * ray.rx.D;

        real3 ryv = ray.ry.P;
        real ty = -(vdot(n, ryv) + d) / vdot(n, ray.ry.D);
        real3 py = ray.ry.P + ty * ray.ry.D;
        
        dPdx = px - p;
        dPdy = py - p;

        // Compute uv deriv.

        int axies[2];

        if (fabs(n[0]) > fabs(n[1]) && fabs(n[0]) > fabs(n[2])) {
            axies[0] = 1; axies[1] = 2;
        } else if (fabs(n[1]) > fabs(n[2])) {
            axies[0] = 0; axies[1] = 2;
        } else {
            axies[0] = 0; axies[1] = 1;
        }
        
        real A[2][2];
        real Bx[2], By[2];
        real x[2], y[2];

        A[0][0] = dPdu[axies[0]];
        A[0][1] = dPdv[axies[0]];
        A[1][0] = dPdu[axies[1]];
        A[1][1] = dPdv[axies[1]];

        Bx[0] = px[axies[0]] - p[axies[0]];
        Bx[1] = px[axies[1]] - p[axies[1]];
        By[0] = py[axies[0]] - p[axies[0]];
        By[1] = py[axies[1]] - p[axies[1]];

        dudx = 1.0; dvdx = 0.0;
        dudy = 0.0; dvdy = 1.0;

        if (SolveLinearSystem22(A, Bx, x)) {
            dudx = x[0]; dvdx = x[1];
        }

        if (SolveLinearSystem22(A, By, y)) {
            dudy = y[0]; dvdy = y[1];
        }

        //printf("p = %f, %f, %f\n", p[0], p[1], p[2]);
        //printf("px = %f, %f, %f\n", px[0], px[1], px[2]);
        //printf("py = %f, %f, %f\n", py[0], py[1], py[2]);
        //printf("dPdx = %f, %f, %f\n", dPdx[0], dPdx[1], dPdx[2]);
        //printf("dPdy = %f, %f, %f\n", dPdy[0], dPdy[1], dPdy[2]);
        //printf("dPdu = %f, %f, %f\n", dPdu[0], dPdu[1], dPdu[2]);
        //printf("dPdv = %f, %f, %f\n", dPdv[0], dPdv[1], dPdv[2]);
        //printf("ddx = %f, %f\n", dudx, dvdx);
        //printf("ddy = %f, %f\n", dudy, dvdy);

    }

    return true;

}

bool ComputeTextureDifferentials(
    real3 &dTdx,                ///< dT/dx [out]
    real3 &dTdy,                ///< dT/dy [out]
    const Ray &ray,
    const real3 &tv0,           ///< Triangle vert
    const real3 &tv1,
    const real3 &tv2,
    const real3 &uv0,           ///< Triangle uv coord
    const real3 &uv1,
    const real3 &uv2)
{

    if (!ray.hasDifferential) {

        dTdx[0] = dTdx[1] = dTdx[2] = 0.0;
        dTdy[0] = dTdy[1] = dTdy[2] = 0.0;

        return false;
    } else {

        //
        // Explicitly find bary coord by casting a ray to the triangle.
        //
        real u0, u1, u2;
        real v0, v1, v2;
        real t0, t1, t2;

        real3 xP, yP;
        real3 xD, yD;

        xP = ray.rx.org;
        yP = ray.ry.org;

        xD = ray.rx.dir;
        yD = ray.ry.dir;

        FindBarycentricCoord(u0, v0, t0, tv0, tv1, tv2, ray.org, ray.dir);
        FindBarycentricCoord(u1, v1, t1, tv0, tv1, tv2, xP, xD);
        FindBarycentricCoord(u2, v2, t2, tv0, tv1, tv2, yP, yD);

        // dTdx = T(px) - T(p) = st1 - st0
        real3 T;
        real3 Tx;
        real3 Ty;

        T  = (1.0 - u0 - v0) * uv0 + u0 * uv1 + v0 * uv2;
        Tx = (1.0 - u1 - v1) * uv0 + u1 * uv1 + v1 * uv2;
        Ty = (1.0 - u2 - v2) * uv0 + u2 * uv1 + v2 * uv2;

        dTdx = Tx - T;
        dTdy = Ty - T;

        //printf("bary:uv0  = %f, %f\n", u0, v0);
        //printf("bary:uv1  = %f, %f\n", u1, v1);
        //printf("bary:uv2  = %f, %f\n", u2, v2);
        //printf("uv0  = %f, %f\n", uv0[0], uv0[1]);
        //printf("uv1  = %f, %f\n", uv1[0], uv1[1]);
        //printf("uv2  = %f, %f\n", uv2[0], uv2[1]);
        //printf("stP  = %f, %f\n", T[0], T[1]);
        //printf("stPx = %f, %f\n", Tx[0], Tx[1]);
        //printf("stPy = %f, %f\n", Ty[0], Ty[1]);
        //printf("dTdx = %f, %f\n", dTdx[0], dTdx[1]);
        //printf("dTdy = %f, %f\n", dTdy[0], dTdy[1]);

    }

    return true;
}
  
