#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>

#include "camera.hh"
#include "kdtree.hh"
#include "parse.hh"
#include "vector.hh"

/*
Vector direct_light(const Ray &r, lights, const Vector &hit, const Material &mat)
{
    Vector light(0, 0, 0);
    float bias = 0.00001f;

    Vector v_normal = r.tri.normal[0];

    Vector shadow_o = r.dir.dot_product(v_normal) <  0 ?
                        hit + bias * v_normal :
                        hit - bias * v_normal; //FIXME normal

    Vector light_amt(0, 0, 0);
    Vector specular_color(0, 0, 0);

    for (const auto &l : lights)
    {
        // Phong model
        Vector dir_light = l.position - r.origin;
        float len = dir_light.dot_product(dir_light);
        dir_light.normalize();

        float l_dot_n = l.dot_product(v_normal);
        if (l_dot_n <= 0)
            continue;

        Ray r_light(shadow_o, dir_light);
        float dist = -1;
        Vector out(0, 0, 0);
        tree.search(r_light, cam, dist, out);
        if (dist == -1 || dist * dist < len)
            continue;

        light_amt *= l.intensity * l_dot_n;

    }

}
*/

int main(int argc, char *argv[])
{
    std::string path_obj = "cube.svati";
    std::string path_mat = "cube.svati";
    std::string path_scene = "cube.svati";

    if (argc > 1)
    {
        path_obj = argv[1];
        if (argc > 2)
            path_mat = argv[2];
        if (argc > 3)
            path_scene = argv[3];
    }

    parse_scene(path_scene);

    int width = 512;
    int height = 512;
    float fov = 120.f;

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
    std::vector<std::string> mat_names;
    mat_names.reserve(map.size());
    for (const auto &it : map)
    {
        //std::cout << it.first << std::endl;
        mat_names.push_back(it.first);
    }

    auto vertices = obj_to_vertices(path_obj, mat_names);
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

    const float ambiant_l[3] = { 0.3f, 0.3f, 0.3f};
    const float dir_l[3] = { 0.2f, 0.2f, 0.2f}; //FIXME

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
            Ray r(o, dir);
            float dist = -1;
            Vector out(0, 0, 0);
            tree.search(r, cam, dist, out);
            if (dist == -1) // not found
                vect[idx] = Vector(0.f, 0.f, 0.f);
            else
            {
                auto material = map.at(mat_names[r.tri.id]);
                /*
                auto direct_l = direct_light(r, d_lights, material);
                Vector indirect_l(0.f, 0.f, 0.f);

                const unsigned nb_ray = 2;
                for (unsigned i = 0; i < nb_ray; ++i)
                {
                    float r1 = static_cast<float>(rand()) / static_cast <float> (RAND_MAX); // fixme random uniform
                    float r2 = static_cast<float>(rand()) / static_cast <float> (RAND_MAX);

                    Vector new_dir(Vector(r1, r2));
                    Ray ray(o, new_dir);

                    tree.search(ray, cam, dist, out, mat); // if intersect, compute indirect_l
                    indirect_l *= 2 * M_PI;

                }
                indirect_l *= static_cast<float>(1 / nb_ray);*/
                vect[idx] = Vector(material.ka.x_ * ambiant_l[0] + material.kd.x_ * dir_l[0] + 0.3,
                                   material.ka.y_ * ambiant_l[1] + material.kd.y_ * dir_l[1],
                                   material.ka.z_ * ambiant_l[2] + material.kd.z_ * dir_l[2]); // Fixme direct light
                /*
                vect[idx] *= (1 / M_PI);
                vect[idx] += indirect_l * 2; // * albedo
                */
            }
        }
    }
    t2 = omp_get_wtime();
    std::cout << "Time raytracing: " << t2 - t1 << "s\n";
    return write_ppm("out.ppm", vect, width, height);

}
