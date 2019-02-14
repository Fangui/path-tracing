#pragma once

#include <cmath>

class Vector
{
public:
    Vector(float x, float y, float z) : x_(x), y_(y), z_(z) {};
    Vector operator+(const Vector &rhs) const;
    Vector operator+=(const Vector &rhs);
    Vector operator*(int lambda) const;
    Vector operator*=(int lambda);
    Vector operator*(const Vector &rhs) const;
    Vector operator*=(const Vector &rhs);
    Vector cross_product(const Vector &rhs) const;
    Vector cross_product_inplace(const Vector &rhs);
    Vector norm(void) const;
    Vector norm_inplace(void);

//private:
    inline float get_dist() { return sqrtf(x_ * x_ + y_ * y_ + z_ * z_); };
    float x_;
    float y_;
    float z_;
};
