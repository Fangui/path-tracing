#include <algorithm>
#include <iostream>
#include <vector>

#include "vector.hh"
#include "camera.hh"
#include "kdtree.hh"

bool compare_y(Vector &a, Vector &b)
{
    return a.y_ < b.y_;
}

int main(void)
{
    /*
    Vector cam_pos(0, 0, -4);
    Vector u(1, 0, 0);
    Vector v(0, 1, 0);

    Camera cam(512, 512, 90, cam_pos, u, v);
    */

    std::vector<Vector> v;

    const int size = 10;
    for (int i = 0; i < size; ++i)
    {
        v.push_back(Vector(i, size - i, (10 - size) % 5));
    }

    KdTree k(v.begin(), v.end(), false);
    /*
    std::sort(v.begin(), v.end(), compare_y);

    for (int i = 0; i < size; ++i)
        std::cout << v[i].x_ << '|';

    std::cout << std::endl;
    */


}
