#pragma once

#include "vector.hh"

//#define EPSILON 1.19209290E-07F
#define EPSILON 0.00001

struct Ray;

struct Triangle
{
    Triangle(Vector &v1, Vector &v2, Vector &v3,
             Vector &n1, Vector &n2, Vector &n3,
             unsigned char id)
    { //FIXME
        vertices[0] = v1;
        vertices[1] = v2;
        vertices[2] = v3;

        normal[0] = n1;
        normal[1] = n2;
        normal[2] = n3;
        this->id = id;

        float x = 0.f;
        float y = 0.f;
        float z = 0.f;

        for (unsigned i = 0; i < 3; ++i)
        {
            x += vertices[i][0];
            y += vertices[i][1];
            z += vertices[i][2];
        }
        mean = Vector(x / 3.f, y / 3.f, z / 3.f);
    }

    Triangle() = default;

    Vector get_mean(void) // return barycentre
    {
        return mean;
    }

    bool intersect(Ray &ray, float &dist) const;
    bool intersect(const Ray &ray) const;

    Vector vertices[3];
    Vector normal[3];
    Vector mean;
       unsigned char id;
};

struct Ray
{
    Ray(Vector &o, Vector &dir) : o(o), dir(dir)
    {
        inv = Vector(1.f / dir[0], 1.f / dir[1], 1.f / dir[2]);
        sign[0] = inv[0] < 0;
        sign[1] = inv[1] < 0;
        sign[2] = inv[2] < 0;
    };
    Vector o;
    Vector dir;
    Vector inv;
    Triangle tri;

    float u;
    float v;

    short sign[3];
};
