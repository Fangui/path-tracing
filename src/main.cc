#include <algorithm>
#include <iostream>
#include <vector>

#include <fstream>

#include "camera.hh"
#include "kdtree.hh"
#include "parse.hh"
#include "vector.hh"

Vector intersect(std::vector<Triangle> &v, const Camera &cam,
                 Vector &dir, Vector &o)
{

    float min_dist = -1;
    int min_idx = -1;

    for (unsigned i = 0; i < v.size(); ++i)
    {
        const Triangle &t = v[i];
        /*
        float dist = -t.normal[0].x * t.vertices[0].x
                    - t.normal[0].y * t.vertices[0].y
                    - t.normal[0].z * t.vertices[0].z;
        */
      //  float tmp = -t[n].dot_product(cam->pos_) + dist;
      //  tmp /= t.normal[0].dot_product(dir);

        Vector out(0, 0, 0);
        if (t.intersect(o, dir, out))
        {
            float distance = (out - cam.pos_).get_dist();
            if (min_dist == -1 || min_dist > distance)
            {
                min_dist = distance;
                min_idx = i;
            }
        }
    }
    if (min_dist != -1)
        return Vector(1, 0, 0);
    return Vector(0, 0, 0);
}

bool inside(const Vector &a,
            const Vector &b,
            const Vector &c,
            const Vector &ray)
{
    Vector b_a = a - b;
    Vector b_c = c - b;
    Vector b_d = ray - b;

    Vector e = b_c.cross_product(b_a);
    Vector f = b_d.cross_product(b_a);

    return e.dot_product(f) > 0;
}

int write_ppm(const std::string &out_path, const std::vector<Vector> &vect,
              int width, int height)
{
    std::ofstream out (out_path);
    unsigned index = 0;
    if (out.is_open())
    {
        out << "P3\n";
        out << width << " " << height << '\n';
        out << 255 << '\n';

        for (int i = 0; i < width; ++i)
        {
            for (int j = 0; j < height; ++j)
            {
                int r = vect[index].x_ * 255.0;
                int g = vect[index].y_ * 255.0;
                int b = vect[index++].z_ * 255.0;
                out << r << " " << g << " " << b << "  ";
            }
            out << '\n';
        }
    }
    else
    {
        std::cerr << "Error while write \n";
        return 1;
    }
    return 0;
}

int main(void)
{
    int width = 512;
    int height = 512;
    float fov = 90.f;

    Vector cam_pos(0, 0, -4);
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

    auto vertices = obj_to_vertices("triangle.svati");
    std::cout << vertices.size() << std::endl;

    std::vector<Vector> vect(width * height);

    Vector C = cam_pos + (w * L); // center

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

            Vector dir = o - cam_pos;

            Triangle tri = vertices[0];
            float dist = -tri.normal[0].x_ * tri.vertices[0].x_
                        - tri.normal[0].y_ * tri.vertices[1].y_
                        - tri.normal[0].z_ * tri.vertices[2].z_;

            float t = -tri.normal[0].dot_product(cam_pos) + dist;
            t /= tri.normal[0].dot_product(dir);

            Vector inter = o * t + cam_pos;


            if (inside(tri.vertices[0], tri.vertices[1], tri.vertices[2], inter)
             && inside(tri.vertices[1], tri.vertices[2], tri.vertices[0], inter)
             && inside(tri.vertices[0], tri.vertices[2], tri.vertices[1], inter))
            {
                float distance = (inter - cam_pos).get_dist();
             //   vect.push_back(Vector(0.8 * 0.65, 0,0));
                vect[idx].set(0.8 * 0.65, 0, 0);
            }
            else
                vect[idx].set(0,0,0);

            //vect.push_back(intersect(vertices, cam, dir, o));
        }
    }
    return write_ppm("out.ppm", vect, width, height);

}
