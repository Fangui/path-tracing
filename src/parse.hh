#pragma once

#include "material.hh"
#include "triangle.hh"

#include <map>

std::vector<Triangle> obj_to_vertices(const std::string &s);
std::map<std::string, Material> parse_materials(const std::string &s);
int write_ppm(const std::string &out_path, const std::vector<Vector> &vect,
              int width, int height);
