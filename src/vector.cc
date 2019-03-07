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

Vector Vector::operator*(float lambda) const
{
    return Vector(tab[0], tab[1], tab[2] ) *= lambda;
}

Vector Vector::operator*=(float lambda)
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

Vector Vector::cross_product(const Vector &rhs) const
{
    return Vector(tab[0], tab[1], tab[2]).cross_product_inplace(rhs);
}

Vector Vector::cross_product_inplace(const Vector &rhs)
{
    float x = tab[1] * rhs[2] - tab[2] * rhs[1];
    float y = tab[2] * rhs[0] - tab[0] * rhs[2];
    float z = tab[0] * rhs[1] - tab[1] * rhs[0];

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
    float dist = this->get_dist();

    for (unsigned i = 0; i < 3; ++i)
        tab[i] /= dist;

    return *this;
}

float Vector::dot_product(const Vector &rhs) const
{
    return tab[0] * rhs[0] + tab[1] * rhs[1] + tab[2] * rhs[2];
}

std::ostream& operator <<(std::ostream& os, const Vector &v)
{
    return os << "x: " << v.tab[0] << " y: " << v.tab[1] << " z: " << v.tab[2];
}
