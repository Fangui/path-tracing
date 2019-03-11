#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
#include <unordered_map>

#include "compute_light.hh"
#include "kdtree.hh"
#include "parse.hh"
#include "vector.hh"

int main(int argc, char *argv[])
{
    std::string path_scene;

    if (argc > 1)
        path_scene = argv[1];
    else
    {
        std::cerr << "Usage: ./main <scene>\n";
        return 1;
    }

    double t1 = omp_get_wtime();

    Scene scene = parse_scene(path_scene);

    Vector u_n = scene.cam_u.norm_inplace();
    Vector v = scene.cam_v.norm_inplace();
    Vector w = v.cross_product(u_n);

    float val = tanf(scene.fov * M_PI / 360);
    val = val == 0.0 ? 0.0001 : val;
    float L = scene.width / 2;
    L /= val; // distance between camera and center of screen

    std::vector<Triangle> vertices;
    for (const auto& name : scene.objs)
      obj_to_vertices(name, scene.mat_names, vertices);

    double t2 = omp_get_wtime();
    std::cout << "Time to parse file: " << t2 - t1 << "s\n";

    t1 = omp_get_wtime();
    auto tree = KdTree(vertices.begin(), vertices.end());
    t2 = omp_get_wtime();
    std::cout << "Time build kdTree: " << t2 - t1 << "s\n";

    std::cout << tree.size() << std::endl;

    std::vector<Vector> vect(scene.width * scene.height);

    Vector C = scene.cam_pos + (w * L); // center

    t1 = omp_get_wtime();

#pragma omp parallel for schedule (dynamic)
    for (int i = -scene.width / 2; i < scene.width / 2; ++i)
    {
        for (int j = -scene.height / 2; j < scene.height / 2; ++j)
        {
            unsigned idx = (i + scene.width / 2) * scene.height + (j + scene.height / 2);
            Vector o = scene.cam_u * j;
            Vector b = scene.cam_v * i;
            o += C;
            o += b;

            Vector dir = o - scene.cam_pos;
            dir.norm_inplace();
            Ray r(scene.cam_pos, dir);

            vect[idx] = cast_ray(scene, r, tree, 0); // depth
        }
    }
    t2 = omp_get_wtime();
    std::cout << "Time raytracing: " << t2 - t1 << "s\n";
    return write_ppm("out.ppm", vect, scene.width, scene.height);

}
