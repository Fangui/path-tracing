#include <algorithm>
#include <iostream>
#include <vector>

#include "camera.hh"
#include "kdtree.hh"
#include "parse.hh"
#include "vector.hh"

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

    auto v = obj_to_vertices("IronMan.obj");

}
