#include "kdtree.hh"
#include "scene.hh"
#include "material.hh"
#include "triangle.hh"

Vector direct_light(const Scene &scene, const Material &mat,
                    const Ray &ray, const KdTree &tree,
                    const Vector& inter, const Vector &normal, int depth);

Vector indirect_light(const Scene &scene,
                      const KdTree &tree, const Vector &inter,
                      const Vector &normal, unsigned char depth);


Vector cast_ray(const Scene &scene, 
                Ray &ray, const KdTree &tree,
                unsigned char depth);
