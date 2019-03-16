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
    std::vector<std::string> mat_names;
    std::unordered_map<std::string, Material> map;
    std::unordered_map<std::string, Texture> map_kd;
};
