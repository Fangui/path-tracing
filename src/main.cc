#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

#include "camera.hh"
#include "kdtree.hh"
#include "parse.hh"
#include "vector.hh"

int main(int argc, char *argv[])
{
    std::string path_obj = "cube.svati";
    std::string path_mat = "cube.svati";

    if (argc > 1)
    {
        path_obj = argv[1];
        if (argc > 2)
            path_mat = argv[2];
    }

    int width = 512;
    int height = 512;
    float fov = 70.f;

    Vector cam_pos(0, 1.5, -5);
    Vector u(1, 0, 0);
    Vector v(0, -1, 0);

    Camera cam(width, height, fov, cam_pos, u, v);

    Vector u_n = u.norm();
    Vector w = v.norm();
    w = w.cross_product_inplace(u_n);

    float val = tanf(cam.fov_ * M_PI / 360);
    val = val == 0.0 ? 0.0001 : val;
    float L = width / 2;
    L /= val; // distance between camera and center of screen

    double t1 = omp_get_wtime();
    auto map = parse_materials(path_mat);
    std::vector<std::string> mat_name;
    mat_name.reserve(map.size());
    for (const auto &it : map)
    {
//        std::cout << it.first << std::endl;
        mat_name.push_back(it.first);
    }

    auto vertices = obj_to_vertices(path_obj, mat_name);
    double t2 = omp_get_wtime();
    std::cout << "Time to parse file: " << t2 - t1 << "s\n";

    t1 = omp_get_wtime();
    auto tree = KdTree(vertices.begin(), vertices.end());
    t2 = omp_get_wtime();
    std::cout << "Time build kdTree: " << t2 - t1 << "s\n";

  //  std::cout << vertices.size() << std::endl;
    std::cout << tree.size() << std::endl;
//    tree.print_infixe();

    std::vector<Vector> vect(width * height);

    Vector C = cam_pos + (w * L); // center

    t1 = omp_get_wtime();
#pragma omp parallel for schedule (dynamic)
    for (int i = -width / 2; i < width / 2; ++i)
    {
        for (int j = -height / 2; j < height / 2; ++j)
        {
            unsigned idx = (i + width / 2) * height + (j + height / 2);
            Vector o = u * j;
            Vector b = v * i;
            o += C;
            o += b;

            Vector dir = cam_pos - o;
            float dist = -1;
            Vector out(0, 0, 0);
            tree.search(o, dir, cam, dist, out);
            if (dist == -1) // not found
                vect[idx] = Vector(0, 0, 0);
            else
                vect[idx] = Vector(1, 0, 0);
        }
    }
    t2 = omp_get_wtime();
    std::cout << "Time raytracing: " << t2 - t1 << "s\n";
    return write_ppm("out.ppm", vect, width, height);

}
