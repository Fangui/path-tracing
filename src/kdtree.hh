#pragma once

#include <functional>
#include <memory>
#include <vector>

#include "vector.hh"

using iterator_v = std::vector<Vector>::iterator;

class KdTree
{
public:
    struct KdNode
    {
        KdNode();
        KdNode(iterator_v beg, iterator_v end);

        std::shared_ptr<KdNode> left;
        std::shared_ptr<KdNode> right;

        float box[6];
    //    Vector &vertex;
        iterator_v beg;
        iterator_v end;
    };

    using childPtr = std::shared_ptr<KdNode>;
    
    KdTree(iterator_v beg, iterator_v end, bool is_vertixes);

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
