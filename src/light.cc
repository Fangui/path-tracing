#include "light.hh"

std::ostream& operator <<(std::ostream& os, const Light &l)
{
    os << "Color :" << l.color << '\n';
    os << "Dir   :" << l.dir;
    return os;
}


Vector Light::compute_light(const Vector &inter, 
                            const KdTree &tree,
                            double &rat) const
{
    Vector l_dir = -dir;

    Vector origin = inter + l_dir * BIAS; // bias
    Ray ray(origin, l_dir);

    if (tree.search_inter(ray))
    {
        rat = 0;
        return l_dir;
    }

    rat = 1;
    return l_dir;
}
