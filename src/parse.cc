#include "kdtree.hh"
#include "parse.hh"
#include "triangle.hh"

#include <fstream>
#include <vector>

inline bool is_blank(char c)
{
    return c == ' ' || c == '\t' || c == '\n';
}

/*
inline bool is_number(char c)
{
    return (c >= '0' && c <= '9') || c == '.';
}*/

std::vector<Triangle> obj_to_vertices(const std::string &s)
{
    std::vector<Triangle> v;
    std::ifstream in(s);

    std::string line;

    float val[9];
    unsigned cpt = 0;
    while (std::getline(in, line))
    {
        if (line[0] == 'v' && line[1] == ' ') // vertices
        {

            for (unsigned i = 2; i < line.size(); ++i)
            {
                std::string s;
                s.reserve(line.size());

                if (is_blank(line[i]))
                    continue;

                while (i < line.size() && !is_blank(line[i]))
                {
                    s += line[i];
                    ++i;
                }
                val[cpt++] = stof(s);
            }
        }
        if (cpt == 9)
        {
            cpt = 0;
            Vector v1(val[0], val[1], val[2]);
            Vector v2(val[3], val[4], val[5]);
            Vector v3(val[6], val[7], val[8]);
            Triangle t(v1, v2, v3);
            v.push_back(t);
        }
    }

    return v;
}
