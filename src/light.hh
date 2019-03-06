#pragma once

#include <iostream>
#include "vector.hh"

struct Light // directional
{
    Light(Vector &color,
          Vector &dir)
    : color(color)
    , dir(dir)
    {}

    Vector color;
    Vector dir;
};

std::ostream& operator <<(std::ostream& os, const Light &l);
