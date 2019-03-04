#pragma once
#include "light.hh"
#include "object.hh"

struct Scene
{
    int height;
    int width;
    Vector cam_pos;
    Vector cam_u;
    Vector cam_v;
    float fov;
    std::vector<std::string> objs;
    std::vector<std::string> mtls;
    std::vector<Light> lights;
    std::vector<Object> objects;
};
