#include "vector.hh"

#include <immintrin.h>

Vector Vector::operator+(const Vector &rhs) const
{
    return Vector(tab[0], tab[1], tab[2] ) += rhs;
}

Vector Vector::operator+=(const Vector &rhs)
{
    tab = _mm256_add_pd(tab, rhs.tab);
    return *this;
}

Vector Vector::operator-(const Vector &rhs) const
{
    return Vector(tab[0], tab[1], tab[2]) -= rhs;
}

Vector Vector::operator-=(const Vector &rhs)
{
    tab = _mm256_sub_pd(tab, rhs.tab);
    return *this;
}

Vector Vector::operator*(double lambda) const
{
    return Vector(tab[0], tab[1], tab[2]) *= lambda;
}

Vector Vector::operator*=(double lambda)
{
    tab = _mm256_mul_pd(tab, _mm256_set1_pd(lambda));
    return *this;
}

Vector Vector::operator*(const Vector &rhs) const
{
    return Vector(tab[0] * rhs[0], tab[1] * rhs[1], tab[2] * rhs[2]);
}

Vector Vector::operator*=(const Vector &rhs)
{
    tab = _mm256_mul_pd(tab, rhs.tab);
    return *this;
}

Vector Vector::operator/=(const Vector &rhs)
{
    tab = _mm256_div_pd(tab, rhs.tab);
    return *this;
}

Vector Vector::operator/(const Vector &rhs) const
{
    return Vector(tab[0] / rhs[0], tab[1] / rhs[1], tab[2] / rhs[2]);
}

Vector Vector::operator/=(double lambda)
{
    tab = _mm256_sub_pd(tab, _mm256_set1_pd(lambda));
    return *this;
}

Vector Vector::operator/(double lambda) const
{
    return Vector(tab[0], tab[1], tab[2]) / lambda;
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

    /*
    tab = _mm256_sub_pd(_mm256_mul_pd(_mm256_shuffle_pd(tab, tab, _MM_SHUFFLE(3, 0, 2, 1)),
                              _mm256_shuffle_pd(rhs.tab, rhs.tab, _MM_SHUFFLE(3,1,0,2))),
                        _mm256_mul_pd(_mm256_shuffle_pd(tab, tab, _MM_SHUFFLE(3, 1, 0, 2)), 
                              _mm256_shuffle_pd(rhs.tab, rhs.tab, _MM_SHUFFLE(3,0,2,1)))); */

    return *this;
}

Vector Vector::norm(void) const
{
    return Vector(tab[0], tab[1], tab[2]).norm_inplace();
}

Vector Vector::norm_inplace(void)
{
    double dist = this->get_dist();
    tab /= dist;

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
