#include <algorithm>
#include <parallel/algorithm>
#include <iostream>
#include <omp.h>

#include "kdtree.hh"

const std::function<bool (Triangle, Triangle)> func[3] =
{
    [](Triangle a, Triangle b) { return a.get_mean().x_ < b.get_mean().x_; },
    [](Triangle a, Triangle b) { return a.get_mean().y_ < b.get_mean().y_; },
    [](Triangle a, Triangle b) { return a.get_mean().z_ < b.get_mean().z_; }
};

#define GET_MIN_MAX(idx, coord) \
    if (box[idx] > beg->vertices[i].coord) \
        box[idx] = beg->vertices[i].coord; \
    if (box[idx + 1] < beg->vertices[i].coord) \
        box[idx + 1] = beg->vertices[i].coord;

static void get_extremum(float box[6], iterator_v beg,
                                       iterator_v end)
{
    box[0] = beg->vertices[0].x_;
    box[1] = beg->vertices[0].x_;

    box[2] = beg->vertices[0].y_;
    box[3] = beg->vertices[0].y_;

    box[4] = beg->vertices[0].z_;
    box[5] = beg->vertices[0].z_;

    ++beg;

    while (beg < end)
    {
        for (unsigned i = 0; i < 3; ++i)
        {
            GET_MIN_MAX(0, x_);
            GET_MIN_MAX(2, y_);
            GET_MIN_MAX(4, z_);
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
    float tmin = (box[ray.sign[0]] - origin.x_) * ray.inv.x_;
    float tmax = (box[1 - ray.sign[0]] - origin.x_) * ray.inv.x_;

    float tymin = (box[2 + ray.sign[1]] - origin.y_) * ray.inv.y_;
    float tymax = (box[3 - ray.sign[1]] - origin.y_) * ray.inv.y_;

    if (tmin > tymax || tymin > tmax)
        return false;

    if (tymin > tmin)
        tmin = tymin;

    if (tymax < tmax)
        tmax = tymax;

    float tzmin = (box[4 + ray.sign[2]] - origin.z_) * ray.inv.z_;
    float tzmax = (box[5 - ray.sign[2]] - origin.z_) * ray.inv.z_;

    if (tmin > tzmax || tzmin > tmax)
        return false;

    return true;
}

#define SEARCH(coord) \
    if (left != nullptr) \
        left.get()->search(ray, cam, dist, last_inter);   \
    if (right != nullptr) \
        right.get()->search(ray, cam, dist, last_inter);

void KdTree::KdNode::search(Ray &ray, const Camera &cam,
                          float &dist, Vector &last_inter)
 {
    if (left == right || inside_box(ray))
    {
        Vector inter(0, 0, 0);
        for (auto it = beg; it < end; ++it)
        {
            if (it->intersect(ray.o, ray.dir, inter))
            {
                float distance = (inter - cam.pos_).get_dist();
                if (dist == -1 || dist > distance)
                {
                    dist = distance;

                    last_inter = inter;
                    ray.tri = *it;
                }
            }
        }

        if (axis == 0)
        {
            SEARCH(x_);
        }
        else if (axis == 1)
        {
            SEARCH(y_);
        }
        else if (axis == 2)
        {
            SEARCH(z_);
        }
    }
}
