#include "vector.hh"

Vector Vector::operator+(const Vector &rhs) const
{
    return Vector(tab[0], tab[1], tab[2] ) += rhs;
}

Vector Vector::operator+=(const Vector &rhs)
{
    for (unsigned i = 0; i < 3; ++i)
        tab[i] += rhs[i];
    return *this;
}

Vector Vector::operator-(const Vector &rhs) const
{
    return Vector(tab[0], tab[1], tab[2]) -= rhs;
}

Vector Vector::operator-=(const Vector &rhs)
{
    for (unsigned i = 0; i < 3; ++i)
        tab[i] -= rhs[i];

    return *this;
}

Vector Vector::operator*(double lambda) const
{
    return Vector(tab[0], tab[1], tab[2]) *= lambda;
}

Vector Vector::operator*=(double lambda)
{
    for (unsigned i = 0; i < 3; ++i)
        tab[i] *= lambda;

    return *this;
}

Vector Vector::operator*(const Vector &rhs) const
{
    return Vector(tab[0] * rhs[0], tab[1] * rhs[1], tab[2] * rhs[2]);
}

Vector Vector::operator*=(const Vector &rhs)
{
    for (unsigned i = 0; i < 3; ++i)
        tab[i] *= rhs[i];

    return *this;
}

Vector Vector::operator/=(const Vector &rhs)
{
    for (unsigned i = 0; i < 3; ++i)
        tab[i] /= rhs[i];

    return *this;
}

Vector Vector::operator/(const Vector &rhs) const
{
    return Vector(tab[0] / rhs[0], tab[1] / rhs[1], tab[2] / rhs[2]);
}

Vector Vector::operator/=(double lambda)
{
    for (unsigned i = 0; i < 3; ++i)
        tab[i] /= lambda;

    return *this;
}

Vector Vector::operator/(double lambda) const
{
    return Vector(tab[0] / lambda, tab[1] / lambda, tab[2] / lambda);
}

Vector operator/(double lambda, const Vector &v)
{
    return v / lambda;
}

Vector operator*(double lambda, const Vector &v)
{
    return v * lambda;
}

Vector Vector::cross_product(const Vector &rhs) const
{
    return Vector(tab[0], tab[1], tab[2]).cross_product_inplace(rhs);
}

Vector Vector::cross_product_inplace(const Vector &rhs)
{
    double x = tab[1] * rhs[2] - tab[2] * rhs[1];
    double y = tab[2] * rhs[0] - tab[0] * rhs[2];
    double z = tab[0] * rhs[1] - tab[1] * rhs[0];

    tab[0] = x;
    tab[1] = y;
    tab[2] = z;

    return *this;
}

Vector Vector::norm(void) const
{
    return Vector(tab[0], tab[1], tab[2]).norm_inplace();
}

Vector Vector::norm_inplace(void)
{
    double dist = this->get_dist();

    for (unsigned i = 0; i < 3; ++i)
        tab[i] /= dist;

    return *this;
}

double Vector::dot_product(const Vector &rhs) const
{
    return tab[0] * rhs[0] + tab[1] * rhs[1] + tab[2] * rhs[2];
}

std::ostream& operator <<(std::ostream& os, const Vector &v)
{
    return os << "x: " << v.tab[0] << " y: " << v.tab[1] << " z: " << v.tab[2];
}
