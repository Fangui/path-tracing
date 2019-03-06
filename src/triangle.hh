#pragma once

#include "vector.hh"

//#define EPSILON 1.19209290E-07F
#define EPSILON 0.00001

struct Triangle
{
    Triangle(Vector &v1, Vector &v2, Vector &v3,
             Vector &n1, Vector &n2, Vector &n3,
             unsigned id)
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
            x += vertices[i].x_;
            y += vertices[i].y_;
            z += vertices[i].z_;
        }
        mean = Vector(x / 3.f, y / 3.f, z / 3.f);
    }

    Triangle() = default;

    Vector get_mean(void) // return barycentre
    {
        return mean;
    }

    bool intersect(const Vector &o,
                   const Vector &ray,
                   Vector& out) const;

    Vector vertices[3];
    Vector normal[3];
    Vector mean;
    unsigned id;
};
