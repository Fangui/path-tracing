#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "camera.hh"
#include "triangle.hh"

using iterator_v = std::vector<Triangle>::iterator;

struct Ray
{
    Ray(Vector o, Vector dir) : o(o), dir(dir)
    {
        inv = Vector(1 / dir.x_, 1 / dir.y_, 1 / dir.z_);
        sign[0] = inv.x_ < 0;
        sign[1] = inv.y_ < 0;
        sign[2] = inv.z_ < 0;
    };
    Vector o;
    Vector dir;
    Vector inv;
    short sign[3];
};

class KdTree
{
public:
    struct KdNode
    {
        KdNode();
        KdNode(iterator_v beg, iterator_v end);

        std::shared_ptr<KdNode> left;
        std::shared_ptr<KdNode> right;

        float box[6]; // pair min : impair max
        iterator_v beg; // beg is median of axis
        iterator_v end; // end = beg + 1 if not leaf
        unsigned char axis = 0;

        void search(Ray &ray, const Camera &cam,
                    float &dist, Vector &last_inter);

        bool inside_box(const Ray &ray) const;

        unsigned size(void)
        {
            unsigned res = std::distance(beg, end);
            if (left)
                res += left.get()->size();
            if (right)
                res += right.get()->size();

            return res;
        }
    };

    using childPtr = std::shared_ptr<KdNode>;

    KdTree(iterator_v beg, iterator_v end);
    void search(const Vector &origin,
                    const Vector &ray, const Camera &cam,
                    float &dist, Vector &last_inter)
    {
        Ray r(origin, ray);
        root_.get()->search(r, cam, dist, last_inter);
    }

    unsigned size(void)
    {
        return root_.get()->size();
    }

private:
    static inline auto make_child()
    {
        return std::make_shared<KdNode>();
    }

    static inline auto make_child(iterator_v beg, iterator_v end)
    {
        return std::make_shared<KdNode>(KdNode(beg, end));
    }

    childPtr root_;
};
