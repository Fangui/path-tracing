#pragma once

#include "vector.hh"

struct Material
{
    Material(double ns, Vector &ka, Vector &kd, Vector &ks,
             Vector &ke, double ni, double d, int illum, Vector &tf,
             const std::string &kd_name,
             const std::string &ka_name)
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
        , ka_name(ka_name)
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
    std::string ka_name;
};
