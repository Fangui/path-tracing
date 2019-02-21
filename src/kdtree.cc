#include <algorithm>
#include <parallel/algorithm>

#include "kdtree.hh"

const std::function<bool (Triangle, Triangle)> func[3] =
{
    [](Triangle a, Triangle b) { return a.get_mean().x_ < b.get_mean().x_; },
    [](Triangle a, Triangle b) { return a.get_mean().y_ < b.get_mean().y_; },
    [](Triangle a, Triangle b) { return a.get_mean().z_ < b.get_mean().z_; }
};

static void get_extremum(float box[6], iterator_v beg,
                                       iterator_v end)
{ //FIXME
    box[0] = (*beg).vertices[0].x_;
    box[1] = (*beg).vertices[0].x_;

    box[2] = (*beg).vertices[0].y_;
    box[3] = (*beg).vertices[0].y_;

    box[4] = (*beg).vertices[0].z_;
    box[5] = (*beg).vertices[0].z_;

    ++beg;
    while (beg < end)
    {
        for (unsigned i = 0; i < 3; ++i) //FIXME
        {
            if (box[0] > beg->vertices[i].x_)
                box[0] = (*beg).vertices[i].x_;

            if (box[1] < (*beg).vertices[i].x_)
                box[1] = (*beg).vertices[i].x_;

            if (box[2] > (*beg).vertices[i].y_)
                box[2] = (*beg).vertices[i].y_;

            if (box[3] < (*beg).vertices[i].y_)
                box[3] = (*beg).vertices[i].y_;

            if (box[4] > (*beg).vertices[i].z_)
                box[4] = (*beg).vertices[i].z_;

            if (box[5] < (*beg).vertices[i].z_)
                box[5] = (*beg).vertices[i].z_;
        }
        ++beg;
    }
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

KdTree::KdTree(iterator_v beg, iterator_v end, bool is_vertices)
{
    root_ = make_child(beg, end);
}

KdTree::KdNode::KdNode(iterator_v beg, iterator_v end)
{
    unsigned dist = std::distance(beg, end);
    if (dist < 20)
    {
        this->beg = beg;
        this->end = end;
        left = nullptr;
        right = nullptr;
    }
    else
    {
        get_extremum(box, beg, end);
        axis = get_longest_axis(box);

        __gnu_parallel::sort(beg, end, func[axis]);
        iterator_v med = beg + dist / 2;

        this->beg = med;
        this->end = med + 1;

        left = make_child(beg, med);
        right = make_child(med, end);
    }
}
