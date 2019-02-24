#include "kdtree.hh"
#include "parse.hh"
#include "triangle.hh"
#include "material.hh"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <map>

inline bool is_blank(char c)
{
    return c == ' ' || c == '\t' || c == '\n';
}

std::map<std::string, Material> parse_materials(const std::string &s)
{
  std::map<std::string, Material> map;
  std::string name;
  float ns, ni, d;
  Vector ka, ks, kd, ke;
  int illum;

  std::ifstream in(s.substr(0, s.length() - 3) + "mtl");
  std::string line;

  while (std::getline(in, line))
  {
    if (line.substr(0, 6) == "newmtl")
    {
      std::string trash;
      name = line.substr(7, line.length());
      in >> trash >> ns;
      in >> trash >> ka.x_ >> ka.y_ >> ka.z_;
      in >> trash >> kd.x_ >> kd.y_ >> kd.z_;
      in >> trash >> ks.x_ >> ks.y_ >> ks.z_;
      in >> trash >> ke.x_ >> ke.y_ >> ke.z_;
      in >> trash >> ni;
      in >> trash >> d;
      in >> trash >> illum;
      Material mat(ns, ka, kd, ks, ke, ni, d, illum);
      // std::cout << "newmtl " << name << std::endl;
      // mat.dump();
      map.emplace(std::make_pair(name, mat));
    }
  }
  return map;
}

std::vector<Triangle> obj_to_vertices(const std::string &s)
{
    parse_materials(s);
    std::vector<Vector> v;
    std::vector<Vector> vn;
    std::vector<Triangle> vt;

    std::ifstream in(s);

    std::string line;

    float val[3];
    while (std::getline(in, line))
    {
        if (line[0] == 'v') // vertices
        {
            unsigned cpt = 0;
            for (unsigned i = 2; i < line.size(); ++i)
            {
                while (i < line.size() && is_blank(line[i]))
                    ++i;

                std::string s;
                s.reserve(line.size() - i);
                while (i < line.size() && !is_blank(line[i]))
                    s += line[i++];

                val[cpt++] = stof(s);
            }

            Vector vect(val[0], val[1], val[2]);
            if (line[1] == 'n') // vn
                vn.push_back(vect);
            else if (line[1] == ' ')
                v.push_back(vect); // v
        }
    }

    vt.reserve(v.size() / 3);
    for (unsigned i = 0; i < v.size() / 3; ++i)
    {
        Triangle t(v[i], v[i + 1], v[i + 2],
                   vn[i], vn[i + 1], vn[i + 2]);

        vt.push_back(t);
    }

    return vt;
}

int write_ppm(const std::string &out_path, const std::vector<Vector> &vect,
              int width, int height)
{
    std::ofstream out (out_path);
    unsigned index = 0;
    if (out.is_open())
    {
        out << "P3\n";
        out << width << " " << height << '\n';
        out << 255 << '\n';

        for (int i = 0; i < width; ++i)
        {
            for (int j = 0; j < height; ++j)
            {
                int r = vect[index].x_ * 255.0;
                int g = vect[index].y_ * 255.0;
                int b = vect[index++].z_ * 255.0;
                out << r << " " << g << " " << b << "  ";
            }
            out << '\n';
        }
    }
    else
    {
        std::cerr << "Error while write in " << out_path << '\n';
        return 1;
    }
    return 0;
}
