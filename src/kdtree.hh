#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "triangle.hh"

using iterator_v = std::vector<Triangle>::iterator;

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

        void search(Ray &ray, const Vector &cam_pos,
                    float &dist);

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
    void search(Ray &r, const Vector &cam_pos,
                    float &dist)
    {
        root_.get()->search(r, cam_pos, dist);
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
