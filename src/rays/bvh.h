
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

    void SAH(const size_t& idx, const size_t& max_leaf_size);
    Trace find_hit(const Ray& ray, const size_t& idx) const;
    Trace hit_queue(const Ray& ray) const;

    void SAH4(const size_t& idx, const size_t& max_leaf_size);
    void bucket_split(const size_t& idx, size_t& rangel, size_t& ranger, BBox& split_leftBox,
                      BBox& split_rightBox);
    Trace hit_queue4(const Ray& ray) const;

    void SAH8(const size_t& idx, const size_t& max_leaf_size);
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

    std::vector<Node4> nodes;

    class Node8 {
        BBox bbox;
        size_t start, size;
        size_t child[8];

        bool is_leaf8() const;
        friend class BVH<Primitive>;
    };
    size_t new_node8(BBox box = {}, size_t start = 0, size_t size = 0);

    // std::vector<Node8> nodes;

    std::vector<Primitive> primitives;
    size_t root_idx = 0;

    class Cluster {
    public:
        BBox representive;
        std::vector<size_t> indexOfPrimitives;

        void updateRepresentive() {
            if(indexOfPrimitives.size()) {
                representive =
                    BBox(m_min /= indexOfPrimitives.size(), m_max /= indexOfPrimitives.size());
            }
        }

        void reset() {
            m_min *= 0;
            m_max *= 0;
        }

        void add(size_t i, BBox bb) {
            indexOfPrimitives.push_back(i);
            m_min += bb.min;
            m_max += bb.max;
        }

    private:
        Vec3 m_min;
        Vec3 m_max;
    };

    struct KBVHNode {
    public:
        BBox bb;
        Cluster c;
        KBVHNode* l;
        KBVHNode* r;

        inline bool isLeaf() const {
            return l == NULL && r == NULL;
        }

        KBVHNode(BBox bb, Cluster c) : bb(bb), c(c), l(NULL), r(NULL) {
        }
    };

    class Kmeans {
    public:
        Kmeans() {
            cluster = NULL;
            children = NULL;
        }
        Kmeans(size_t iterCount, size_t K, size_t P, std::vector<size_t> indexOfPrimitives) {
            m_iterations = iterCount;
            m_K = K;
            m_P = P;
            this->indexOfPrimitives = indexOfPrimitives;
            cluster = new Cluster[K];
            children = new Kmeans*[K];
            for(size_t i = 0; i < K; i++) {
                children[i] = NULL;
            }
        };
        size_t m_iterations;
        size_t m_K;
        size_t m_P;

        Cluster* cluster;
        std::vector<size_t> indexOfPrimitives;

        Kmeans** children;

        KBVHNode* root;
    };
    void constructKaryTree(Kmeans* kmeans, size_t max_leaf_size);
    void allocateKaryTree(Kmeans* kmeans, size_t max_leaf_size);
    Kmeans* kmeans_root;

    Trace hit_kmeans(const Ray& ray, Kmeans* kmeans) const;
};

} // namespace PT

#ifdef CARDINAL3D_BUILD_REF
#include "../reference/bvh.inl"
#else
#include "../student/bvh.inl"
#endif
