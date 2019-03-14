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

        double box[6]; // pair min : impair max
        iterator_v beg; // beg is median of axis
        iterator_v end; // end = beg + 1 if not leaf
        unsigned char axis = 0;

        void search(Ray &ray, double &dist) const;

        bool search_inter(const Ray &ray) const;

        bool inside_box(const Ray &ray) const;

        inline bool is_child(void) const { return left == right; }

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

            std::cout << "extremum: ";
            for (unsigned i = 0; i < 6; ++i)
                std::cout << box[i] << " ";
            std::cout << '\n';
            for (auto it = beg; it < end; ++it)
            {
                std::cout << it->get_mean() << '\n';
            }
            if (right)
                right.get()->print_infixe();

            /*
            for (auto it = beg; it < end; ++it)
            {
                for (int i = 0; i < 3; ++i)
                {
                    if (box[0] >= it[0].vertices[i][0])
                    {
                        std::cout << box[0] << std::endl;
                        std::cout << it[0].vertices[i][0] << std::endl;
                        exit(20);
                    }
                    if (box[1] <= it[0].vertices[i][0])
                    {
                        std::cout << "min" << box[0] << std::endl;
                        std::cout << "max:"<<  box[1] << std::endl;
                        std::cout << it[0].vertices[i][0] << std::endl;
                            
                        std::cout << std::distance(beg, end) << std::endl;
                        for (auto it2 = beg; it2 < end; ++it2)
                        {
                            for (int j = 0; j < 3; ++j)
                            {
                                std::cout << it2->vertices[j][0] << std::endl;
                            }
                        }
                       exit(30);
                    }
                    if (box[2] >= it[0].vertices[i][1])
                        exit(4);

                    if (box[3] <= it[0].vertices[i][1])
                        exit(5);

                    if (box[4] >= it[0].vertices[i][2])
                        exit(6);

                    if (box[5] <= it[0].vertices[i][2])
                    {
                        std::cout << box[5] << std::endl;
                        std::cout << it[0].vertices[i][2] << std::endl;
                        exit(7);
                    }
                }
            }*/

            if (left)
                left.get()->print_infixe();

            if (right)
                right.get()->print_infixe();

        }
    };

    using childPtr = std::unique_ptr<KdNode>;

    KdTree(iterator_v beg, iterator_v end);
    bool search(Ray &r, double &dist) const
{
        root_.get()->search(r, dist);

        return dist != -1;
    }


    bool search_inter(const Ray &r) const
    {
        return root_.get()->search_inter(r);
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
