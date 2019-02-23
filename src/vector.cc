#include "vector.hh"

Vector Vector::operator+(const Vector &rhs) const
{
    return Vector(x_, y_, z_ ) += rhs;
}

Vector Vector::operator+=(const Vector &rhs)
{
    x_ += rhs.x_;
    y_ += rhs.y_;
    z_ += rhs.z_;

    return *this;
}

Vector Vector::operator-(const Vector &rhs) const
{
    return Vector(x_, y_, z_ ) -= rhs;
}

Vector Vector::operator-=(const Vector &rhs)
{
    x_ -= rhs.x_;
    y_ -= rhs.y_;
    z_ -= rhs.z_;

    return *this;
}

Vector Vector::operator*(float lambda) const
{
    return Vector(x_, y_, z_ ) *= lambda;
}

Vector Vector::operator*=(float lambda)
{
    x_ *= lambda;
    y_ *= lambda;
    z_ *= lambda;

    return *this;
}

Vector Vector::operator*(const Vector &rhs) const
{
    return Vector(x_ * rhs.x_, y_ * rhs.y_, z_ * rhs.z_);
}

Vector Vector::operator*=(const Vector &rhs)
{
    x_ *= rhs.x_;
    y_ *= rhs.y_;
    z_ *= rhs.z_;

    return *this;
}

Vector Vector::cross_product(const Vector &rhs) const
{
    return Vector(x_, y_, z_).cross_product_inplace(rhs);
}

Vector Vector::cross_product_inplace(const Vector &rhs)
{
    float x = y_ * rhs.z_ - z_ * rhs.y_;
    float y = z_ * rhs.x_ - x_ * rhs.z_;
    float z = x_ * rhs.y_ - y_ * rhs.x_;

    x_ = x;
    y_ = y;
    z_ = z;

    return *this;
}

Vector Vector::norm(void) const
{
    return Vector(x_, y_, z_).norm_inplace();
}

Vector Vector::norm_inplace(void)
{
    float dist = this->get_dist();

    x_ /= dist;
    y_ /= dist;
    z_ /= dist;

    return *this;
}

float Vector::dot_product(const Vector &rhs) const
{
    return x_ * rhs.x_ + y_ * rhs.y_ + z_ * rhs.z_;
}

std::ostream& operator <<(std::ostream& os, const Vector &v)
{
    return os << "x: " << v.x_ << " y: " << v.y_ << " z: " << v.z_;
}
