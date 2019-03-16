#pragma once

#include <iostream>
#include "kdtree.hh"
#include "triangle.hh"
#include "vector.hh"

struct Light // directional
{
    Light(const Vector &color,
          const Vector &dir)
    : color(color)
    , dir(dir)
    {}

    virtual ~Light() = default;

    virtual Vector compute_light(const Vector &inter,
                          const KdTree &tree,
                          double &rat) const;

    Vector color;
    Vector dir;
};

std::ostream& operator <<(std::ostream& os, const Light &l);
