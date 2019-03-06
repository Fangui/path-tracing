#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>
#include <omp.h>
#include <unordered_map>

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

    Camera cam(scene.width, scene.height, scene.fov, scene.cam_pos, scene.cam_u, scene.cam_v);

    Vector u_n = scene.cam_u.norm();
    Vector w = scene.cam_v.norm();
    w = w.cross_product_inplace(u_n);

    float val = tanf(scene.fov * M_PI / 360);
    val = val == 0.0 ? 0.0001 : val;
    float L = scene.width / 2;
    L /= val; // distance between camera and center of screen

    std::unordered_map<std::string, Material> map;
    for (const auto& name : scene.mtls)
        parse_materials(name, map);
    std::vector<std::string> mat_names;
    mat_names.reserve(map.size());
    for (const auto &it : map)
    {
        //std::cout << it.first << std::endl;
        mat_names.push_back(it.first);
    }

    std::vector<Triangle> vertices;
    for (const auto& name : scene.objs)
      obj_to_vertices(name, mat_names, vertices);
    double t2 = omp_get_wtime();
    std::cout << "Time to parse file: " << t2 - t1 << "s\n";

    t1 = omp_get_wtime();
    auto tree = KdTree(vertices.begin(), vertices.end());
    t2 = omp_get_wtime();
    std::cout << "Time build kdTree: " << t2 - t1 << "s\n";

  //  std::cout << vertices.size() << std::endl;
    std::cout << tree.size() << std::endl;
//    tree.print_infixe();

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

            Vector dir = scene.cam_pos - o;
            Ray r(o, dir);
            float dist = -1;
            Vector out(0, 0, 0);
            tree.search(r, cam, dist, out);
            if (dist == -1) // not found
                vect[idx] = Vector(0.f, 0.f, 0.f);
            else
            {
                auto material = map.at(mat_names[r.tri.id]);
                vect[idx] = Vector(material.ka.x_ * scene.a_light.x_,
                                   material.ka.y_ * scene.a_light.y_,
                                   material.ka.z_ * scene.a_light.z_); // ambient light
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
                
                /*
                vect[idx] *= (1 / M_PI);
                vect[idx] += indirect_l * 2; // * albedo
                */
            }
        }
    }
    t2 = omp_get_wtime();
    std::cout << "Time raytracing: " << t2 - t1 << "s\n";
    return write_ppm("out.ppm", vect, scene.width, scene.height);

}
