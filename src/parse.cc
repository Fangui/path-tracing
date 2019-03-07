#include "kdtree.hh"
#include "parse.hh"
#include "triangle.hh"
#include "light.hh"
#include "material.hh"
#include "json.hpp"

#include <iostream>
#include <algorithm>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <exception>
#include <string>

using json = nlohmann::json;

inline bool is_separator(char c)
{
    return c == ' ' || c == '\t' || c == '\n' || c == '/';
}

Scene parse_scene(const std::string& filename)
{
    Scene scene;
    std::ifstream in(filename);
    if (!in)
        std::cerr << "Json parser error: file " + filename + " not found";
    json j;
    in >> j;

    try{
        scene.height = (j["height"].get<int>());
        scene.width = (j["width"].get<int>());

        auto cpos = j["camera"]["pos"].begin();
        scene.cam_pos.set(cpos[0], cpos[1], cpos[2]);

        cpos = j["camera"]["u"].begin();
        scene.cam_u.set(cpos[0], cpos[1], cpos[2]);

        scene.fov = j["camera"]["fov"];

        cpos = j["camera"]["v"].begin();
        scene.cam_v.set(cpos[0], cpos[1], cpos[2]);

        auto objs = j["meshs"];
        for (auto e : objs)
            scene.objs.emplace_back(e);

        auto mtls = j["materials"];
        for (auto e : mtls)
            scene.mtls.emplace_back(e);

        auto lights = j["lights"];
        for (auto e : lights)
        {
            std::string s = e["type"].get<std::string>();
            int r = 1;
            Vector color(1, 1, 1);
            Vector dir(1, 1, 1);

            if (e.find("color") != e.end())
            {
                auto pos = e["color"].begin();
                color.set(pos[0], pos[1], pos[2]);
            }

            if (e.find("vector") != e.end())
            {
                auto pos = e["vector"].begin();
                dir.set(pos[0], pos[1], pos[2]);
            }


            if (s == "ambient")
                scene.a_light = color;
            else if (s == "directional")
            {
                scene.lights.push_back(Light(color, dir));
            }
            else
            {
                std::cerr << "Light not implemented yet" << std::endl;
                exit(2);
            }
        }

        auto objects = j["objects"];
        for (auto e : objects)
        {
            int mesh = e["mesh"].get<int>();
            int mtl = e["mtl"].get<int>();

            auto pos = e["position"].begin();
            Vector p(1, 1, 1);
            p.set(pos[0], pos[1], pos[2]);

            auto rot = e["rotation"].begin();
            Vector r(1, 1, 1);
            r.set(rot[0], rot[1], rot[2]);

            auto sc = e["scale"].begin();
            Vector s(1, 1, 1);
            s.set(sc[0], sc[1], sc[2]);
            scene.objects.emplace_back(Object(mesh, mtl, p, r, s));
        }

        //Dump scene
        std::cout << "----- Scene -----" << std::endl;
        std::cout << "Height : " << scene.height << std::endl;
        std::cout << "Width : " << scene.width << std::endl;
        std::cout << "Camera :" << std::endl;
        std::cout << "     pos : " << scene.cam_pos << '\n';
        std::cout << "     fov : " << scene.fov << std::endl;
        std::cout  << std::endl;
        std::cout << "Meshs :" << std::endl;
        for (auto e : scene.objs)
            std::cout << "     " << e << std::endl;
        std::cout  << std::endl;
        std::cout << "Materials :" << std::endl;
        for (auto e : scene.mtls)
            std::cout << "     " << e << std::endl;
        std::cout  << std::endl;

        std::cout << "Lights :" << std::endl;
        std::cout << scene.a_light << "\n\n";
        std::cout << "Directional_light" << '\n';
        for (auto e : scene.lights)
            std::cout<< e << '\n';
        std::cout  << std::endl;
        std::cout << "Objects :" << std::endl;
        for (auto e : scene.objects)
        {
            std::cout << "     mesh/mtl : " << e.mesh << " " << e.mtl << std::endl;
            std::cout << "     pos : " << e.pos[0] << " " << e.pos[1] << " " << e.pos[2] << std::endl;
            std::cout << "     rot : " << e.rot[0] << " " << e.rot[1] << " " << e.rot[2] << std::endl;
            std::cout << "     scale : " << e.scale[0] << " " << e.scale[1] << " " << e.scale[2] << std::endl;
            std::cout  << std::endl;
        }
    } catch (std::exception& e){
        std::cout << e.what() << std::endl;
    }
    return scene;
}
void parse_materials(const std::string &s, std::unordered_map<std::string, Material>& map)
{
    std::string name;
    std::ifstream in(s);
    std::string line;

    if (!in.is_open())
    {
        std::cerr << "Materiel not found: " << s << '\n';
        exit(2);
    }

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
                    strin >> trash >> ka.get_x() >> ka.get_y() >> ka.get_z();
                else if (id == "Kd")
                    strin >> trash >> kd.get_x() >> kd.get_y() >> kd.get_z();
                else if (id == "Ks")
                    strin >> trash >> ks.get_x() >> ks.get_y() >> ks.get_z();
                else if (id == "Ke")
                    strin >> trash >> ke.get_x() >> ke.get_y() >> ke.get_z();
                else if (id == "Ni")
                    strin >> trash >> ni;
                else if (id == "d ")
                    strin >> trash >> d;
                else if (id == "il")
                    strin >> trash >> illum;
                else if (id == "Tf")
                    strin >> trash >> tf.get_x() >> tf.get_y() >> tf.get_z();
                else if (id == "ma")
                    continue;
                else
                    break;
                id.clear();
                id += in.peek();
                if (id == "n")
                    break;
            }
            Material mat(ns, ka, kd, ks, ke, ni, d, illum, tf);
            //  std::cout << "newmtl " << name << std::endl;
            //  mat.dump();
            map.emplace(std::make_pair(name, mat));
        }
    }
}

void obj_to_vertices(const std::string &s, const std::vector<std::string> &mat_names,
                     std::vector<Triangle>& v_tri)
{
    std::vector<Vector> v;
    std::vector<Vector> vn;

    std::ifstream in(s);

    std::string line;
    std::unordered_map<std::string, unsigned> map;

    float val[3];
    unsigned idx[9] = { 0 };
    unsigned cur_idx = 0;

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
            auto name = line.substr(7, line.length());

            unsigned cpt = 0;
            for (const auto &str : mat_names)
            {
                if (str == name)
                    break;
                ++cpt;
            }


            if (cpt == mat_names.size())
                std::cerr << "Material name " << name << " not found \n";

            cur_idx = cpt;
        }
        else if (line[0] == 'f')
        {
            unsigned cpt = 0;
            for (unsigned i = 2; i < line.size(); ++i)
            {
                if (i < line.size() && line[i] == '/')
                {
                    i += 1;
                    ++cpt;
                }

                while (i < line.size() && is_separator(line[i]))
                    ++i;

                std::string s;
                s.reserve(line.size() - i);
                while (i < line.size() && !is_separator(line[i]))
                    s += line[i++];

                idx[cpt++] = stof(s) - 1;
            } // FIXME
            if (vn.size() == 0) // if not vn in file
            {
                Triangle t(v[idx[0]], v[idx[3]], v[idx[6]],
                           v[0], v[0], v[0], cur_idx);
                v_tri.push_back(t);
            }
            else
            {
                Triangle t(v[idx[0]], v[idx[3]], v[idx[6]],
                           vn[idx[2]], vn[idx[5]], vn[idx[8]], cur_idx);
                v_tri.push_back(t);
            }
        }
    }
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
                int r = vect[index][0] * 255.0;
                int g = vect[index][1] * 255.0;
                int b = vect[index++][2] * 255.0;
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
