#pragma once

#include <cmath>
#include <iostream>

class Vector
{
public:
    Vector() = default;
    Vector(float x, float y, float z) : x_(x), y_(y), z_(z) {};
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

//private:
    float x_;
    float y_;
    float z_;
};

std::ostream& operator <<(std::ostream& os, const Vector &v);
