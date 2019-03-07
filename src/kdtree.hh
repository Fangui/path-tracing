#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "camera.hh"
#include "triangle.hh"

using iterator_v = std::vector<Triangle>::iterator;

struct Ray
{
    Ray(Vector &o, Vector &dir) : o(o), dir(dir)
    {
        inv = Vector(1.f / dir[0], 1.f / dir[1], 1.f / dir[2]);
        sign[0] = inv[0] < 0;
        sign[1] = inv[1] < 0;
        sign[2] = inv[2] < 0;


    };
    Vector o;
    Vector dir;
    Vector inv;
    Triangle tri;
    short sign[3];
};

class KdTree
{
public:
    struct KdNode
    {
        KdNode();
        KdNode(iterator_v beg, iterator_v end);

        std::unique_ptr<KdNode> left;
        std::unique_ptr<KdNode> right;

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

        void print_infixe(void)
        {
            if (left)
                left.get()->print_infixe();

            std::cout << "etremum: ";
            for (unsigned i = 0; i < 6; ++i)
                std::cout << box[i] << " ";
            std::cout << '\n';
            for (auto it = beg; it < end; ++it)
            {
                std::cout << it->get_mean() << '\n';
            }
            if (right)
                right.get()->print_infixe();
        }
    };

    using childPtr = std::unique_ptr<KdNode>;

    KdTree(iterator_v beg, iterator_v end);
    void search(Ray &r, const Camera &cam,
                    float &dist, Vector &last_inter)
    {
        root_.get()->search(r, cam, dist, last_inter);
    }

    void print_infixe()
    {
        root_.get()->print_infixe();
    }

    unsigned size(void)
    {
        return root_.get()->size();
    }

private:
    static inline auto make_child()
    {
        return std::make_unique<KdNode>();
    }

    static inline auto make_child(iterator_v beg, iterator_v end)
    {
        return std::make_unique<KdNode>(KdNode(beg, end));
    }

    childPtr root_;
};
