#include "kdtree.hh"
#include "scene.hh"
#include "material.hh"
#include "triangle.hh"

Vector direct_light(const Material &mat, const Scene &scene,
                    const Ray &ray, const KdTree &tree,
                    const Vector& out);
