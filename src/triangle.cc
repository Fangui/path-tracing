#include "triangle.hh"

bool Triangle::intersect(const Vector &o,
                         const Vector &ray,
                         Vector &out) const
{
    const float EPSILON = 0.0000001;

    Vector vertex0 = vertices[0];
    Vector vertex1 = vertices[1];
    Vector vertex2 = vertices[2];

    Vector edge1 = vertex1 - vertex0;
    Vector edge2 = vertex2 - vertex0;
    Vector h = ray.cross_product(edge2);

    float det = edge1.dot_product(h); //FIXME
    if (det > -EPSILON && det < EPSILON)
        return false;    // This ray is parallel to this triangle.
    float f = 1.f / det;
    Vector s = o - vertex0;
    float u = f * (s.dot_product(h));

    if (u < 0.0 || u > 1.0)
        return false;

    Vector q = s.cross_product(edge1);
    float v = f * (ray.dot_product(q));
    if (v < 0.0 || u + v > 1.0)
        return false;

    // At this stage we can compute t to find out where the intersection point is on the line.
    float t = f * edge2.dot_product(q);
    if (t > EPSILON) // ray intersection
    {
        out = o + ray * t;
        return true;
    }
    return false;
}
