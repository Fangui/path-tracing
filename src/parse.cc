#include "kdtree.hh"
#include "parse.hh"
#include "triangle.hh"

#include <iostream>
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
    std::vector<Vector> v;
    std::vector<Vector> vn;
    std::vector<Triangle> vt;

    std::ifstream in(s);

    std::string line;

    float val[3];
    unsigned cpt = 0;
    while (std::getline(in, line))
    {
        cpt = 0;
        if (line[0] == 'v') // vertices
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

        if (cpt == 3)
        {
            Vector vect(val[0], val[1], val[2]);
            if (line[1] == 'n') // vn
                vn.push_back(vect);
            else if (line[1] == ' ')
                v.push_back(vect); // v
        }
    }

    vt.reserve(v.size() / 3);
    for (unsigned i = 0; i < v.size(); i += 3)
    {
        unsigned idx = i / 3;
        Triangle t(v[idx], v[idx + 1], v[idx + 2], 
                   vn[idx], vn[idx + 1], vn[idx + 2]);

        vt.push_back(t);
    }

    return vt;
}
