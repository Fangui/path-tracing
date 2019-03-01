#pragma once

struct Light
{
    Light(const std::string& s, Vector p, int r)
    : type(s)
    , pos(p)
    , radius(r)
    {}

    std::string type;
    Vector pos;
    int radius;
};
