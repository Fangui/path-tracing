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

void denoise(std::vector<Vector>& vect, std::vector<Vector>& res,
             unsigned s, unsigned width)
{
    std::vector<float> vectr;
    std::vector<float> vectg;
    std::vector<float> vectb;
    unsigned len = vect.size();

    for (unsigned x = width; x < len - width; ++x)
    {
        vectr.clear();
        vectb.clear();
        vectg.clear();
        for (unsigned i = 0; i < s; i++)
        {
            for (unsigned j = 0; j < s; j++)
            {
                unsigned ind = (i - s / 2 + 1) * width + (j - s / 2 + 1);
                if (x + ind < len)
                {
                    vectr.emplace_back(vect[x + ind][0]);
                    vectg.emplace_back(vect[x + ind][1]);
                    vectb.emplace_back(vect[x + ind][2]);
                }
            }
        }
        std::sort(vectr.begin(), vectr.end());
        std::sort(vectg.begin(), vectg.end());
        std::sort(vectb.begin(), vectb.end());
        res.emplace_back(Vector(vectr[s * s / 2 - 1],
                                vectg[s * s / 2 - 1], vectb[s * s / 2 - 1]));
    }
}

int main(int argc, char *argv[])
{
    std::string path_scene;
    std::string out_file = "out";
    int matrix_size = 3;

    if (argc > 1)
        path_scene = argv[1];
    else
    {
        std::cerr << "Usage: ./main <scene> <nb_ray> <depth> <filter> <out_file>\n";
        return 1;
    }


    double t1 = omp_get_wtime();

    Scene scene = parse_scene(path_scene);
    if (argc > 2)
    {
        scene.nb_ray = atoi(argv[2]);
        if (argc > 3)
            scene.depth = atoi(argv[3]);

        if (argc > 4)
            matrix_size = atoi(argv[4]);

        if (argc > 5)
            out_file = argv[5];
    }

    Vector u_n = scene.cam_u.norm_inplace();
    Vector v = scene.cam_v.norm_inplace();
    Vector w = v.cross_product(u_n);

    double val = tanf(scene.fov * M_PI / 360);
    val = val == 0.0 ? 0.0001 : val;
    double L = scene.width / 2;
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

//    constexpr double gamma = 1. / 2.2;
#pragma omp parallel for schedule (dynamic)
    for (int i = -scene.width / 2; i < scene.width / 2; ++i)
    {
        for (int j = scene.height / 2; j > -scene.height / 2; --j)
        {
            unsigned idx = (i + scene.width / 2) * scene.height + (scene.height / 2 - j);
            Vector o = scene.cam_u * j;
            Vector b = scene.cam_v * i;
            o += C;
            o += b;

            Vector dir = o - scene.cam_pos;
            dir.norm_inplace();
            Ray r(scene.cam_pos, dir);

            vect[idx] = cast_ray(scene, r, tree, 0); // depth
//            for (unsigned g = 0; g < 3; ++g) // gamme
//                vect[idx][g] = pow(vect[idx][g], gamma);
        }
    }
    t2 = omp_get_wtime();
    std::cout << "Time raytracing: " << t2 - t1 << "s\n";

    write_ppm(out_file + ".ppm", vect, scene.width, scene.height);
    double t4 = omp_get_wtime();
    /*
    std::vector<Vector> out;

    for (int i = 0; i < scene.width; i += 2)
    {
        for (int j = 0; j < scene.height; j += 2)
        {
            Vector c = (vect[i * scene.width + j]
                      + vect[i * scene.width + j + 1]
                      + vect[(i + 1) * scene.height + j]
                      + vect[(i + 1) * scene.height + j + 1]) / 4;

            out.push_back(c);
        }
    }*/
    std::vector<Vector> res;
    denoise(vect, res, matrix_size, scene.width);

    double t3 = omp_get_wtime();
    std::cout << "Time applying denoise: " << t3 - t4 << "s\n";

    //return write_ppm("out.ppm", out, scene.width / 2, scene.height / 2);
    return write_ppm(out_file + "_denoise.ppm", res, scene.width, scene.height);

}
