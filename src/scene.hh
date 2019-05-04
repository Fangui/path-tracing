#pragma once

#include "light.hh"
#include "material.hh"
#include "matrix.hh"
#include "object.hh"
#include "sphere_light.hh"
#include "texture.hh"

#include <vector>
#include <unordered_map>

struct Scene
{
    ~Scene()
    {
        for (unsigned i = 0; i < lights.size(); ++i)
            delete lights[i];
    }

    int nb_ray = 16;
    unsigned char depth = 2;
    int height;
    int width;
    Vector cam_pos;
    Vector cam_u;
    Vector cam_v;
    double fov;
    std::vector<std::string> objs;
    std::vector<std::string> mtls;
    Vector a_light;
    Matrix transform;
    std::vector<Light*> lights;
    std::vector<Object> objects;
    std::vector<Vector> emissive; // coord of emmisive lights
    std::vector<std::string> emissive_name;

    std::vector<std::string> mat_names;
    std::unordered_map<std::string, Material> map;
    std::unordered_map<std::string, Texture> map_text;
};
