#pragma once

#include "light.hh"

struct SphereLight : public Light
{
    SphereLight(const Vector &color,
                const Vector &dir,
                int r)
        : Light(color, dir), 
          r(r)
    {}

    ~SphereLight() = default;

    Vector get_pos() const
    {
        return dir;
    }

    Vector compute_light(const Vector &inter, 
                         const KdTree &tree,
                         double &rat) const override;

    int r;
};
