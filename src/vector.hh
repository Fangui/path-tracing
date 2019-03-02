#pragma once

#include <cmath>
#include <iostream>

class Vector
{
public:
    Vector() : x_(0), y_(0), z_(0) {};
    Vector(float x, float y, float z) : x_(x), y_(y), z_(z) {};
    Vector(float r1, float r2) // create a hemisphere
    {
        float sinTheta = sqrtf(1 - r1 * r1);
        float phi = 2 * M_PI * r2;
        float x = sinTheta * cosf(phi);
        float z = sinTheta * sinf(phi);
        x_ = x;
        y_ = phi;
        z_ = z;
    };
    Vector operator+(const Vector &rhs) const;
    Vector operator+=(const Vector &rhs);
    Vector operator-(const Vector &rhs) const;
    Vector operator-=(const Vector &rhs);
    Vector operator*(float lambda) const;
    Vector operator*=(float lambda);
    Vector operator*(const Vector &rhs) const;
    Vector operator*=(const Vector &rhs);
    Vector cross_product(const Vector &rhs) const;
    Vector cross_product_inplace(const Vector &rhs);
    Vector norm(void) const;
    Vector norm_inplace(void);

    float dot_product(const Vector &rhs) const;

    float get_dist() { return sqrtf(x_ * x_ + y_ * y_ + z_ * z_); };

    void set(float x, float y, float z)
    {
        x_ = x;
        y_ = y;
        z_ = z;
    }

    Vector reflect(const Vector &i, const Vector &n)
    {
        return i - n * i.dot_product(n) * 2; // check formule
 //       return I - 2 * dotProduct(I, N) * N;
    }
//private:
    float x_;
    float y_;
    float z_;
};

std::ostream& operator <<(std::ostream& os, const Vector &v);
