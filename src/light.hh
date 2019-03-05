#pragma once

struct Light
{
    Light(Vector &color,
          Vector &dir)
    : color(color)
    , dir(dir)
    {}

    Vector color;
    Vector dir;
};
