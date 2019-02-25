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
        id = 0;
        sign[0] = inv.x_ < 0;
        sign[1] = inv.y_ < 0;
        sign[2] = inv.z_ < 0;
    };
    Vector o;
    Vector dir;
    Vector inv;
    unsigned id;
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

    using childPtr = std::shared_ptr<KdNode>;

    KdTree(iterator_v beg, iterator_v end, std::vector<std::string> &mat_names);
    void search(Ray &r, const Camera &cam,
                    float &dist, Vector &last_inter,
                    std::string &mat)
    {
        root_.get()->search(r, cam, dist, last_inter);
        if (dist != -1)
            mat = mat_names[r.id];
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
        return std::make_shared<KdNode>();
    }

    static inline auto make_child(iterator_v beg, iterator_v end)
    {
        return std::make_shared<KdNode>(KdNode(beg, end));
    }

    childPtr root_;
    std::vector<std::string> mat_names;
};
