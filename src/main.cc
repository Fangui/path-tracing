#include "vector.hh"
#include "camera.hh"

int main(void)
{
    Vector cam_pos(0, 0, -4);
    Vector u(1, 0, 0);
    Vector v(0, 1, 0);

    Camera cam(512, 512, 90, cam_pos, u, v);
}
