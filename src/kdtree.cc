#include <algorithm>
#include <parallel/algorithm>
#include <iostream>
#include <omp.h>

#include "kdtree.hh"

const std::function<bool (Triangle, Triangle)> func[3] =
{
    [](Triangle a, Triangle b) { return a.get_mean()[0] < b.get_mean()[0]; },
    [](Triangle a, Triangle b) { return a.get_mean()[1] < b.get_mean()[1]; },
    [](Triangle a, Triangle b) { return a.get_mean()[2] < b.get_mean()[2]; }
};

#define GET_MIN_MAX(idx, coord) \
    if (box[idx] > beg->vertices[i][coord]) \
        box[idx] = beg->vertices[i][coord]; \
    if (box[idx + 1] < beg->vertices[i][coord]) \
        box[idx + 1] = beg->vertices[i][coord];

static void get_extremum(float box[6], iterator_v beg,
                                       iterator_v end)
{
    box[0] = beg->vertices[0][0];
    box[1] = beg->vertices[0][0];

    box[2] = beg->vertices[0][1];
    box[3] = beg->vertices[0][1];

    box[4] = beg->vertices[0][2];
    box[5] = beg->vertices[0][2];

    ++beg;

    while (beg < end)
    {
        for (unsigned i = 0; i < 3; ++i)
        {
            GET_MIN_MAX(0, 0);
            GET_MIN_MAX(2, 1);
            GET_MIN_MAX(4, 2);
        }
        ++beg;
    }

    for (unsigned i = 0; i < 6; i += 2) // expand bounding box
        box[i] -= 0.1;
    for (unsigned i = 1; i < 6; i += 2)
        box[i] += 0.1;
}

static unsigned get_longest_axis(float box[6])
{
    float diff_x = fabs(box[1] - box[0]);
    float diff_y = fabs(box[3] - box[2]);
    float diff_z = fabs(box[5] - box[4]);

    if (diff_x > diff_y)
    {
        if (diff_x > diff_z)
            return 0;
        return 2;
    }
    if (diff_y > diff_z)
        return 1;
    return 2;
}

KdTree::KdTree(iterator_v beg, iterator_v end)
{
    root_ = make_child(beg, end);
}

KdTree::KdNode::KdNode(iterator_v beg, iterator_v end)
{
    unsigned dist = std::distance(beg, end);
    get_extremum(box, beg, end);
    if (dist < 4)
    {
        this->beg = beg;
        this->end = end;
        left = nullptr;
        right = nullptr;
    }
    else
    {
        axis = get_longest_axis(box);

        __gnu_parallel::sort(beg, end, func[axis]);
        iterator_v med = beg + dist / 2;

        left = make_child(beg, med);
        right = make_child(med + 1, end);

        this->beg = med;
        this->end = med + 1;
    }
}

bool KdTree::KdNode::inside_box(const Ray &ray) const
{
    const Vector &origin = ray.o;
    float tmin = (box[ray.sign[0]] - origin[0]) * ray.inv[0];
    float tmax = (box[1 - ray.sign[0]] - origin[0]) * ray.inv[0];

    float tymin = (box[2 + ray.sign[1]] - origin[1]) * ray.inv[1];
    float tymax = (box[3 - ray.sign[1]] - origin[1]) * ray.inv[1];

    if (tmin > tymax || tymin > tmax)
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (box[4 + ray.sign[2]] - origin[2]) * ray.inv[2];
    float tzmax = (box[5 - ray.sign[2]] - origin[2]) * ray.inv[2];

    if (tmin > tzmax || tzmin > tmax)
        return false;

    return true;
}

void KdTree::KdNode::search(Ray &ray, const Vector &cam_pos,
                          float &dist, Vector &last_inter)
 {
    if (left == right || inside_box(ray))
    {
        Vector inter(0, 0, 0);
        for (auto it = beg; it < end; ++it)
        {
            if (it->intersect(ray.o, ray.dir, inter))
            {
                float distance = fabs(inter[2] - cam_pos[2]);
                if (dist > distance || dist == -1)
                {
                    dist = distance;

                    last_inter = inter;
                    ray.tri = *it;
                }
            }
        }

        if (axis == 2) // Good opti but may be dangerous
        {
            if (cam_pos[2] < beg->vertices[2][2])
            {
                float tmp = dist;
                left.get()->search(ray, cam_pos, dist, last_inter);
                if (tmp != dist)
                    return;
                right.get()->search(ray, cam_pos, dist, last_inter);
            }
            else
            {
                float tmp = dist;
                right.get()->search(ray, cam_pos, dist, last_inter);
                if (tmp != dist)
                    return;
                left.get()->search(ray, cam_pos, dist, last_inter);

            }
        }
        else
        {
            if (left != nullptr)
                left.get()->search(ray, cam_pos, dist, last_inter);

            if (right != nullptr)
                right.get()->search(ray, cam_pos, dist, last_inter);
        }
    }
}
