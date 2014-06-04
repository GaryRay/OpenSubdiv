#ifndef __MALLIE_CAMERA_H__
#define __MALLIE_CAMERA_H__

#include <string>
#include "common.h"

class Camera {
public:
    void BuildCameraFrame(double const eye[3],
                          double const lookat[3],
                          double const up[3],
                          double fov,
                          int width, int height);

    Ray GenerateRay(double u, double v) const;

    double eye_[3];
    double up_[3];
    double lookat_[3];

    // In world space
    double origin_[3];
    double corner_[3];
    double du_[3];
    double dv_[3];
    double fov_;
};

#endif // __LIINA_CAMERA_H__

// vim:set sw=2 ts=2 expandtab:
