#pragma once

#include <cmath>
#include <iostream>

inline float to_rad(int deg)
{
    return deg * (M_PI / 180.0);
}

class Vector
{
public:
    Vector()
    {
        for (unsigned i = 0; i < 3; ++i)
            tab[i] = 0;
        tab[3] = 1;
    };
    Vector(float x, float y, float z)
    {
      tab[0] = x;
      tab[1] = y;
      tab[2] = z;
      tab[3] = 1;
    }

    /*
    Vector(float r1, float r2) // create a hemisphere
    {
        float sinTheta = sqrtf(1 - r1 * r1);
        float phi = 2 * M_PI * r2;
        float x = sinTheta * cosf(phi);
        float z = sinTheta * sinf(phi);
        x_ = x;
        y_ = phi;
        z_ = z;
    };*/

    Vector operator+(const Vector &rhs) const;
    Vector operator+=(const Vector &rhs);
    Vector operator-(const Vector &rhs) const;
    Vector operator-=(const Vector &rhs);
    Vector operator*(float lambda) const;
    Vector operator*=(float lambda);
    Vector operator*(const Vector &rhs) const;
    Vector operator*=(const Vector &rhs);
    Vector operator/(const Vector &rhs) const;
    Vector operator/=(const Vector &rhs);
    Vector operator/(float lambda) const;
    Vector operator/=(float lambda);

    float operator[](unsigned idx) const { return tab[idx]; };
    float& operator[](unsigned idx) { return tab[idx]; };

   // Vector operator*(float lambda, const Vector &rhs);

    Vector cross_product(const Vector &rhs) const;
    Vector cross_product_inplace(const Vector &rhs);
    Vector norm(void) const;
    Vector norm_inplace(void);

    float dot_product(const Vector &rhs) const;

    float get_dist() { return sqrtf(tab[0] * tab[0] + tab[1] * tab[1] + tab[2] * tab[2]); };

    void set(float x, float y, float z)
    {
        tab[0] = x;
        tab[1] = y;
        tab[2] = z;
    }

    friend std::ostream& operator <<(std::ostream& os, const Vector &v);

private:
    float tab[4];
};

Vector operator/(float lambda, const Vector &v);
Vector operator*(float lambda, const Vector &v);
