#ifndef __MALLIE_RAYDIFF_H__
#define __MALLIE_RAYDIFF_H__

#include "common.h"

namespace mallie {

///
/// Compute differentials for a transfer ray against screen coordinate (x, y).
/// @return true If derivatives was successfully computed.
/// @return false If derivatives was failed to compute.
///
bool ComputeDifferentialsForTransfer(
        real3 &dPdx,                ///< dP/dx [out]
        real3 &dPdy,                ///< dP/dy [out]
        real      &dudx,                ///< du/dx [out]
        real      &dvdx,                ///< dv/dx [out]
        real      &dudy,                ///< du/dy [out]
        real      &dvdy,                ///< dv/dy [out]
        const Ray &ray,
        const real3 &p,
        const real3 &n,
        const real3 &dPdu,          ///< dP/du (tangent)
        const real3 &dPdv);         ///< dP/dv (binormal)

///
/// Compute differentials for a reflection ray against screen coordinate (x, y).
/// @return true If derivatives was successfully computed.
/// @return false If derivatives was failed to compute.
///
bool ComputeDifferentialsForReflection(
    real3       &dDdx,                ///< dD/dx [out]               
    real3       &dDdy,                ///< dD/dy [out]
    const Ray       &ray,
    const real3 &P,
    const real3 &D,
    const real3 &N,
    const real3 &p0,
    const real3 &p1,
    const real3 &p2,
    const real3 &n0,
    const real3 &n1,
    const real3 &n2);

    ///
    /// Compute differentials for a refraction ray against screen coordinate (x, y).
    /// @return true If derivatives was successfully computed.
    /// @return false If derivatives was failed to compute.
    ///
    bool ComputeDifferentialsForRefraction(
        real3       &dDdx,                  ///< dD/dx [out]
        real3       &dDdy,                  ///< dD/dy [out]
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
        real             eta);

    ///
    /// Compute texture coord differentials against screen coordinate (x, y).
    /// @return true If derivatives was successfully computed.
    /// @return false If derivatives was failed to compute.
    ///
    bool ComputeTextureDifferentials(
        real3 &dTdx,                ///< dT/dx [out]
        real3 &dTdy,                ///< dT/dy [out]
        const Ray &ray,
        const real3 &v0,            ///< Triangle vert
        const real3 &v1,
        const real3 &v2,
        const real3 &uv0,           ///< Triangle uv coord
        const real3 &uv1,
        const real3 &uv2);

}; // namespace


#endif // __MALLIE_RAYDIFF_H__
