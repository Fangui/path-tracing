#pragma once

#include <cmath>
#include <iostream>

inline double to_rad(int deg)
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
    Vector(double x, double y, double z)
    {
      tab[0] = x;
      tab[1] = y;
      tab[2] = z;
      tab[3] = 1;
    }

    /*
    Vector(double r1, double r2) // create a hemisphere
    {
        double sinTheta = sqrtf(1 - r1 * r1);
        double phi = 2 * M_PI * r2;
        double x = sinTheta * cosf(phi);
        double z = sinTheta * sinf(phi);
        x_ = x;
        y_ = phi;
        z_ = z;
    };*/

    Vector operator+(const Vector &rhs) const;
    Vector operator+=(const Vector &rhs);
    Vector operator-(const Vector &rhs) const;
    Vector operator-=(const Vector &rhs);
    Vector operator*(double lambda) const;
    Vector operator*=(double lambda);
    Vector operator*(const Vector &rhs) const;
    Vector operator*=(const Vector &rhs);
    Vector operator/(const Vector &rhs) const;
    Vector operator/=(const Vector &rhs);
    Vector operator/(double lambda) const;
    Vector operator/=(double lambda);

    double operator[](unsigned idx) const { return tab[idx]; };
    double& operator[](unsigned idx) { return tab[idx]; };

   // Vector operator*(double lambda, const Vector &rhs);

    Vector cross_product(const Vector &rhs) const;
    Vector cross_product_inplace(const Vector &rhs);
    Vector norm(void) const;
    Vector norm_inplace(void);

    double dot_product(const Vector &rhs) const;

    double get_dist() { return sqrtf(tab[0] * tab[0] + tab[1] * tab[1] + tab[2] * tab[2]); };

    void set(double x, double y, double z)
    {
        tab[0] = x;
        tab[1] = y;
        tab[2] = z;
    }

    friend std::ostream& operator <<(std::ostream& os, const Vector &v);

private:
    double tab[4];
};

Vector operator/(double lambda, const Vector &v);
Vector operator*(double lambda, const Vector &v);
