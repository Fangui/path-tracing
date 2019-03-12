#pragma once

#include "vector.hh"

struct Material
{
    Material(double ns, Vector &ka, Vector &kd, Vector &ks,
             Vector &ke, double ni, double d, int illum, Vector &tf,
             const std::string &kd_name)
        : ns(ns)
        , ka(ka)
        , kd(kd)
        , ks(ks)
        , ke(ke)
        , ni(ni)
        , d(d)
        , illum(illum)
        , tf(tf)
        , kd_name(kd_name)
    { }

    void dump();

    double ns;
    Vector ka;
    Vector kd;
    Vector ks;
    Vector ke;
    double ni;
    double d;
    int illum;
    Vector tf;
    std::string kd_name;
};
