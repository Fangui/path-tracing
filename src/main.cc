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
        return Vector(255, 0, 0);
    return Vector(0, 0, 0);
}

int main(void)
{
    int width = 512;
    int height = 512;

    Vector cam_pos(0, 0, -4);
    Vector u(1, 0, 0);
    Vector v(0, 1, 0);

    Camera cam(width, height, 90, cam_pos, u, v);

    u = u.norm_inplace();
    v = v.norm_inplace();

    Vector w = v.cross_product(u);

    float val = tanf(cam.fov_ * M_PI / 360);
    val = val == 0.0 ? 0.0001 : val;
    float L = width / 2;
    L /= val;

    auto vertices = obj_to_vertices("triangle.svati");

    std::vector<Vector> vect;
    vect.reserve(width * height);

    Vector C = cam_pos + (w * L);

    for (int i = -width / 2; i < width / 2; ++i)
    {
        for (int j = -height / 2; j < height / 2; ++j)
        {

            Vector a = u * i;
            Vector b = v * j;

            Vector o = a + C + b;
            Vector dir = o - cam_pos;

            vect.push_back(intersect(vertices, cam, dir, o));
        }
    }

    std::ofstream out ("out.ppm");
    unsigned index = 0;
    if (out.is_open())
    {
        out << "P3\n";
        out << width << " " << height << '\n';
        out << 255 << '\n';

        for (int i = -width / 2; i < width / 2; ++i)
        {
            for (int j = -height / 2; j < height / 2; ++j)
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
        std::cerr << "Error while write \n";
    return 0;
}
