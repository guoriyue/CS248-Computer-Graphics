
#pragma once

#include "../lib/mathlib.h"
#include "../platform/gl.h"

#include "trace.h"

namespace PT {

template<typename Primitive> class BVH {
public:
    BVH() = default;
    BVH(std::vector<Primitive>&& primitives, size_t max_leaf_size = 1);
    void build(std::vector<Primitive>&& primitives, size_t max_leaf_size = 1);

    BVH(BVH&& src) = default;
    BVH& operator=(BVH&& src) = default;

    BVH(const BVH& src) = delete;
    BVH& operator=(const BVH& src) = delete;

    BBox bbox() const;
    Trace hit(const Ray& ray) const;

    BVH copy() const;
    size_t visualize(GL::Lines& lines, GL::Lines& active, size_t level, const Mat4& trans) const;

    std::vector<Primitive> destructure();
    void clear();

    void SAH(size_t idx, size_t max_leaf_size);
    Trace find_hit(const Ray& ray, size_t idx) const;
    Trace hit_queue(const Ray& ray) const;
    
    void SAH4(size_t idx, size_t max_leaf_size);
    void bucket_split(size_t idx, size_t& rangel, size_t& ranger, BBox& split_leftBox, BBox& split_rightBox);
    Trace hit_queue4(const Ray& ray) const;

    void SAH8(size_t idx, size_t max_leaf_size);
    Trace hit_queue8(const Ray& ray) const;

private:
    class Node {
        BBox bbox;
        size_t start, size, l, r;

        bool is_leaf() const;
        friend class BVH<Primitive>;
    };
    size_t new_node(BBox box = {}, size_t start = 0, size_t size = 0, size_t l = 0, size_t r = 0);

    // std::vector<Node> nodes;

    class Node4 {
        BBox bbox;
        size_t start, size;
        size_t child[4];

        bool is_leaf4() const;
        friend class BVH<Primitive>;
    };
    size_t new_node4(BBox box = {}, size_t start = 0, size_t size = 0);

    // std::vector<Node4> nodes;


    class Node8 {
        BBox bbox;
        size_t start, size;
        size_t child[8];

        bool is_leaf8() const;
        friend class BVH<Primitive>;
    };
    size_t new_node8(BBox box = {}, size_t start = 0, size_t size = 0);

    std::vector<Node8> nodes;

    std::vector<Primitive> primitives;
    size_t root_idx = 0;
};

} // namespace PT

#ifdef CARDINAL3D_BUILD_REF
#include "../reference/bvh.inl"
#else
#include "../student/bvh.inl"
#endif
