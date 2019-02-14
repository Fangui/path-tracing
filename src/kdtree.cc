#include <algorithm>
#include <parallel/algorithm>

#include "kdtree.hh"

const std::function<bool (Vector, Vector)> func[3] = 
{
    [](Vector a, Vector b) { return a.x_ < b.x_; },
    [](Vector a, Vector b) { return a.y_ < b.y_; },
    [](Vector a, Vector b) { return a.z_ < b.z_; }
};

static void get_extremum(float box[6], iterator_v beg,
                                iterator_v end)
{ //FIXME
    box[0] = (*beg).x_;
    box[1] = (*beg).x_;

    box[2] = (*beg).y_;
    box[3] = (*beg).y_;

    box[4] = (*beg).z_;
    box[5] = (*beg).z_;

    ++beg;
    while (beg < end)
    {
        if (box[0] > (*beg).x_)
            box[0] = (*beg).x_;

        if (box[1] < (*beg).x_)
            box[1] = (*beg).x_;

        if (box[2] > (*beg).y_)
            box[2] = (*beg).y_;

        if (box[3] < (*beg).y_)
            box[3] = (*beg).y_;

        if (box[4] > (*beg).z_)
            box[4] = (*beg).z_;

        if (box[5] < (*beg).z_)
            box[5] = (*beg).z_;
        ++beg;
    }
}

static void get_mean_vertices(std::vector<Vector> &v, 
                       iterator_v beg,
                       iterator_v end)
{
    while (beg < end)
    {
        float x = ((*beg).x_ + (*(beg + 1)).x_ + (*(beg + 2)).x_) / 3.f;
        float y = ((*beg).y_ + (*(beg + 1)).y_ + (*(beg + 2)).y_) / 3.f;
        float z = ((*beg).z_ + (*(beg + 1)).z_ + (*(beg + 2)).z_) / 3.f;
        v.push_back(Vector(x, y, z));
        beg += 3;
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
    if (is_vertices)
    {
        unsigned dist = std::distance(beg, end);
        std::vector<Vector> v;
        v.reserve(dist / 3);

        get_mean_vertices(v, beg, end);
    }
    else
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
        unsigned axis = get_longest_axis(box);

        __gnu_parallel::sort(beg, end, func[axis]);
        iterator_v med = beg + dist / 2;

        this->beg = med;
        this->end = med + 1;

        left = make_child(beg, med);
        right = make_child(med, end);
    }
}
