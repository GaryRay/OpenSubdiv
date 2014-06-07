#define _USE_MATH_DEFINES
#include <cmath>
#include <cstdio>
#include <algorithm>

#include "camera.h"

static inline double vdot(double const a[3], double const b[3]) {
  return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

static inline void vcross(double c[3], double const a[3], double const b[3]) {
  c[0] = a[1] * b[2] - a[2] * b[1];
  c[1] = a[2] * b[0] - a[0] * b[2];
  c[2] = a[0] * b[1] - a[1] * b[0];
}

static inline double vlength(double const v[3]) {
  double len2 = vdot(v, v);
  if (std::abs(len2) > 1.0e-30) {
    return sqrt(len2);
  }
  return 0.0;
}

static void vnormalize(double v[3]) {
  float len = vlength(v);
  if (std::abs(len) > 1.0e-30) {
    double inv_len = 1.0 / len;
    v[0] *= inv_len;
    v[1] *= inv_len;
    v[2] *= inv_len;
  }
}

void Camera::BuildCameraFrame(double const eye[3], double const lookat[3],
                              double const up[3], double fov,
                              int width, int height)
{
    double flen =
        (0.5f * (double)height / tanf(0.5f * (double)(fov * M_PI / 180.0f)));
    double look1[3];
    look1[0] = lookat[0] - eye[0];
    look1[1] = lookat[1] - eye[1];
    look1[2] = lookat[2] - eye[2];
    // vcross(u, up1, look1);
    // flip
    double u[3], v[3], corner[3], origin[3];
    vcross(u, look1, up);
    vnormalize(u);

    vcross(v, look1, u);
    vnormalize(v);

    double aspect = height/(float)width;
    u[0] *= aspect;
    u[1] *= aspect;
    u[2] *= aspect;

    v[0] /= aspect;
    v[1] /= aspect;
    v[2] /= aspect;

    vnormalize(look1);
    look1[0] = flen * look1[0] + eye[0];
    look1[1] = flen * look1[1] + eye[1];
    look1[2] = flen * look1[2] + eye[2];
    corner[0] = look1[0] - 0.5f * (width * u[0] + height * v[0]);
    corner[1] = look1[1] - 0.5f * (width * u[1] + height * v[1]);
    corner[2] = look1[2] - 0.5f * (width * u[2] + height * v[2]);

    origin[0] = eye[0];
    origin[1] = eye[1];
    origin[2] = eye[2];

    // Store intermediate
    origin_[0] = origin[0];
    origin_[1] = origin[1];
    origin_[2] = origin[2];

    corner_[0] = corner[0];
    corner_[1] = corner[1];
    corner_[2] = corner[2];

    du_[0] = u[0];
    du_[1] = u[1];
    du_[2] = u[2];

    dv_[0] = v[0];
    dv_[1] = v[1];
    dv_[2] = v[2];

    fov_ = fov;
}

Ray Camera::GenerateRay(double u, double v) const {
  real3 dir;

  dir[0] = (corner_[0] + u * du_[0] + v * dv_[0]) - origin_[0];
  dir[1] = (corner_[1] + u * du_[1] + v * dv_[1]) - origin_[1];
  dir[2] = (corner_[2] + u * du_[2] + v * dv_[2]) - origin_[2];
  dir.normalize();

  real3 org;
  org[0] = origin_[0];
  org[1] = origin_[1];
  org[2] = origin_[2];

  Ray ray;
  ray.org = org;
  ray.dir = dir;

  return ray;
}
