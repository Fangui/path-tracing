#pragma once

#include "vector.hh"

struct Triangle
{
    Triangle(Vector &a, Vector &b, Vector &c)
    { //FIXME
        vertices[0] = a;
        vertices[1] = b;
        vertices[2] = c;
    }

    /*: vertices[0](a),
                                vertices[1](b), vertices[2](c) {};
*/

    Vector get_mean() //FIXME
    {
        float x = 0.f;
        float y = 0.f;
        float z = 0.f;

        for (unsigned i = 0; i < 3; ++i)
        {
            x += vertices[i].x_;
            y += vertices[i].y_;
            z += vertices[i].z_;
        }

        return Vector(x / 3.f, y / 3.f, z / 3.f);
    }

    Vector vertices[3];
};
