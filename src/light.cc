#include "light.hh"

std::ostream& operator <<(std::ostream& os, const Light &l)
{
    os << "Color :" << l.color << '\n';
    os << "dir :" << l.dir;
    return os;
}
