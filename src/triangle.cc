#include "triangle.hh"

bool Triangle::intersect(Ray &ray,
                         double &dist) const
{
    const Vector &vertex0 = vertices[0];
    const Vector &vertex1 = vertices[1];
    const Vector &vertex2 = vertices[2];

    const Vector edge1 = vertex1 - vertex0;
    const Vector edge2 = vertex2 - vertex0;
    Vector h = ray.dir.cross_product(edge2);

    double det = edge1.dot_product(h);
    if (det > -EPSILON && det < EPSILON)
        return false;    // This ray is parallel to this triangle.
    double f = 1.f / det;
    Vector s = ray.o - vertex0;
    double u = f * (s.dot_product(h));

    if (u < 0.0 || u > 1.0)
        return false;

    s.cross_product_inplace(edge1);
    double v = f * (ray.dir.dot_product(s));
    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot_product(s);
    if (t > EPSILON) // ray intersection
    {
        ray.u = u;
        ray.v = v;
        dist = t;
        return true;
    }
    return false;
}

bool Triangle::intersect(const Ray &ray) const
{
    Vector vertex0 = vertices[0];
    Vector vertex1 = vertices[1];
    Vector vertex2 = vertices[2];

    Vector edge1 = vertex1 - vertex0;
    Vector edge2 = vertex2 - vertex0;
    Vector h = ray.dir.cross_product(edge2);

    double det = edge1.dot_product(h);
    if (det > -EPSILON && det < EPSILON)
        return false;    // This ray is parallel to this triangle.
    double f = 1.f / det;
    Vector s = ray.o - vertex0;
    double u = f * (s.dot_product(h));

    if (u < 0.0 || u > 1.0)
        return false;

    s = s.cross_product_inplace(edge1);
    double v = f * (ray.dir.dot_product(s));
    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    double t = f * edge2.dot_product(s);
    if (t > EPSILON) // ray intersection
        return true;
    
    return false;
}
