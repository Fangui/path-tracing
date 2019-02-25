#include "kdtree.hh"
#include "parse.hh"
#include "triangle.hh"
#include "material.hh"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>

inline bool is_separator(char c)
{
    return c == ' ' || c == '\t' || c == '\n' || c == '/';
}

std::map<std::string, Material> parse_materials(const std::string &s)
{
  std::map<std::string, Material> map;
  std::string name;
//  std::ifstream in(s.substr(0, s.length() - 3) + "mtl");
  std::ifstream in(s);
  std::string line;

  while (std::getline(in, line))
  {
    float ns = 0;
    float ni = 0;
    float d = 0;
    int illum = 0;
    Vector ka, ks, kd, ke, tf;
    if (line.substr(0, 6) == "newmtl")
    {
      std::string trash;
      name = line.substr(7, line.length());
      while (std::getline(in, line))
      {
        auto id = line.substr(0, 2);
        std::istringstream strin(line);
        if (id == "Ns")
          strin >> trash >> ns;
        else if (id == "Ka")
          strin >> trash >> ka.x_ >> ka.y_ >> ka.z_;
        else if (id == "Kd")
          strin >> trash >> kd.x_ >> kd.y_ >> kd.z_;
        else if (id == "Ks")
          strin >> trash >> ks.x_ >> ks.y_ >> ks.z_;
        else if (id == "Ke")
          strin >> trash >> ke.x_ >> ke.y_ >> ke.z_;
        else if (id == "Ni")
          strin >> trash >> ni;
        else if (id == "d ")
          strin >> trash >> d;
        else if (id == "il")
          strin >> trash >> illum;
        else if (id == "Tf")
          strin >> trash >> tf.x_ >> tf.y_ >> tf.z_;
        else if (id == "ma")
          continue;
        else
          break;
        id.clear();
        id += in.peek();
        if (id == "n")
          break;
      }
/*      in >> trash >> ns;
      in >> trash >> ka.x_ >> ka.y_ >> ka.z_;
      in >> trash >> kd.x_ >> kd.y_ >> kd.z_;
      in >> trash >> ks.x_ >> ks.y_ >> ks.z_;
      in >> trash >> ke.x_ >> ke.y_ >> ke.z_;
      in >> trash >> ni;
      in >> trash >> d;
      in >> trash >> illum;*/
      Material mat(ns, ka, kd, ks, ke, ni, d, illum, tf);
    //  std::cout << "newmtl " << name << std::endl;
    //  mat.dump();
      map.emplace(std::make_pair(name, mat));
    }
  }
  return map;
}

std::vector<Triangle> obj_to_vertices(const std::string &s,
                                      const std::vector<std::string> &mat_names)
{
    std::vector<Vector> v;
    std::vector<Vector> vn;
    std::vector<Triangle> v_tri;

    std::ifstream in(s);

    std::string line;
    std::map<std::string, unsigned> map;

    float val[3];
    unsigned idx[9] = { 0 };
    unsigned idx_cur;

    while (std::getline(in, line))
    {
        if (line[0] == 'v') // vertices
        {
            unsigned cpt = 0;
            for (unsigned i = 2; i < line.size(); ++i)
            {
                while (i < line.size() && is_separator(line[i]))
                    ++i;

                std::string s;
                s.reserve(line.size() - i);
                while (i < line.size() && !is_separator(line[i]))
                    s += line[i++];

                val[cpt++] = stof(s);
            }

            Vector vect(val[0], val[1], val[2]);
            if (line[1] == 'n') // vn
                vn.push_back(vect);
            else if (line[1] == ' ')
                v.push_back(vect); // v
        }
        else if (line.substr(0, 6) == "usemtl")
        {
            const auto name = line.substr(7, line.length());
            unsigned cpt = 0;
            for (const auto &str : mat_names)
            {
                if (str == name)
                    break;
                ++cpt;
            }

            if (cpt == mat_names.size())
                std::cerr << "Material name not found \n";
            idx_cur = cpt;
        }
        else if (line[0] == 'f')
        {
            unsigned cpt = 0;
            for (unsigned i = 2; i < line.size(); ++i)
            {
                while (i < line.size() && is_separator(line[i]))
                    ++i;

                std::string s;
                s.reserve(line.size() - i);
                while (i < line.size() && !is_separator(line[i]))
                    s += line[i++];

                idx[cpt++] = stof(s);
            } // FIXME
            Triangle t(v[idx[0]], v[idx[3]], v[idx[6]],
                    vn[idx[2]], vn[idx[5]], vn[idx[8]]);

            v_tri.push_back(t);
        }
}

//v_tri.reserve(v.size() / 3);


    return v_tri;
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
