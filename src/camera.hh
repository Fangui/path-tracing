#pragma once

#include "vector.hh"

struct Camera
{
    Camera(int width, int height, float fov,
           Vector &pos, Vector &u,
           Vector &v) : width_(width), height_(height),
                                   fov_(fov), pos_(pos), u_(u), v_(v) {};

    int width_;
    int height_;
    float fov_;
    Vector &pos_;
    Vector &u_;
    Vector &v_;
};
