#include "kdtree.hh"
#include "parse.hh"
#include "triangle.hh"
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
        scene.cam_pos.x_ = *(cpos++);
        scene.cam_pos.y_ = *(cpos++);
        scene.cam_pos.z_ = *(cpos++);

        cpos = j["camera"]["u"].begin();
        scene.cam_u.x_ = *(cpos++);
        scene.cam_u.y_ = *(cpos++);
        scene.cam_u.z_ = *(cpos++);
        scene.fov = j["camera"]["fov"];

        cpos = j["camera"]["v"].begin();
        scene.cam_v.x_ = *(cpos++);
        scene.cam_v.y_ = *(cpos++);
        scene.cam_v.z_ = *(cpos++);
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
            Vector p(1, 1, 1);
            if (e.find("vector") != e.end())
            {
                auto pos = e["vector"].begin();
                p.x_ = *(pos++);
                p.y_ = *(pos++);
                p.z_ = *(pos++);
            }
            scene.lights.emplace_back(Light(s, p, r));
        }

        auto objects = j["objects"];
        for (auto e : objects)
        {
            int mesh = e["mesh"].get<int>();
            int mtl = e["mtl"].get<int>();

            auto pos = e["position"].begin();
            Vector p(1, 1, 1);
            p.x_ = *(pos++);
            p.y_ = *(pos++);
            p.z_ = *(pos++);

            auto rot = e["rotation"].begin();
            Vector r(1, 1, 1);
            r.x_ = *(rot++);
            r.y_ = *(rot++);
            r.z_ = *(rot++);

            auto sc = e["scale"].begin();
            Vector s(1, 1, 1);
            s.x_ = *(sc++);
            s.y_ = *(sc++);
            s.z_ = *(sc++);
            scene.objects.emplace_back(Object(mesh, mtl, p, r, s));
        }

        //Dump scene
        std::cout << "----- Scene -----" << std::endl;
        std::cout << "Height : " << scene.height << std::endl;
        std::cout << "Width : " << scene.width << std::endl;
        std::cout << "Camera :" << std::endl;
        std::cout << "     pos : [ " << scene.cam_pos.x_ << ", " <<  scene.cam_pos.y_ << ", " <<
            scene.cam_pos.z_ << " ]" << std::endl;
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
        for (auto e : scene.lights)
            std::cout << "     " << e.type << " " << e.pos.x_ <<
                " " << e.pos.y_ <<
                " " << e.pos.z_ << std::endl;
        std::cout  << std::endl;
        std::cout << "Objects :" << std::endl;
        for (auto e : scene.objects)
        {
            std::cout << "     mesh/mtl : " << e.mesh << " " << e.mtl << std::endl;
            std::cout << "     pos : " << e.pos.x_ << " " << e.pos.y_ << " " << e.pos.z_ << std::endl;
            std::cout << "     rot : " << e.rot.x_ << " " << e.rot.y_ << " " << e.rot.z_ << std::endl;
            std::cout << "     scale : " << e.scale.x_ << " " << e.scale.y_ << " " << e.scale.z_ << std::endl;
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

            auto e = name.find(':');
            if (e != std::string::npos)
            {
                name = name.substr(e + 1, name.length());
                std::cout << name << std::endl;
            }
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
                while (i < line.size() && is_separator(line[i]))
                    ++i;

                std::string s;
                s.reserve(line.size() - i);
                while (i < line.size() && !is_separator(line[i]))
                    s += line[i++];

                idx[cpt++] = stof(s) - 1;
            } // FIXME
            Triangle t(v[idx[0]], v[idx[3]], v[idx[6]],
                    v[idx[2]], v[idx[5]], v[idx[8]], cur_idx); // FIXME replace v by vn

            v_tri.push_back(t);
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
