#pragma once

#include "vector.hh"

struct Material
{
    Material(float ns, Vector &ka, Vector &kd, Vector &ks,
             Vector &ke, float ni, float d, int illum)
        : ns(ns)
        , ka(ka)
        , kd(kd)
        , ks(ks)
        , ke(ke)
        , ni(ni)
        , d(d)
        , illum(illum)
    { }

    void dump();

    float ns;
    Vector ka;
    Vector kd;
    Vector ks;
    Vector ke;
    float ni;
    float d;
    int illum;
};
