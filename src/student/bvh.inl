
#include "../rays/bvh.h"
#include "debug.h"
#include <stack>

#include "../rays/ispc_bvh.h"

namespace PT {

// construct BVH hierarchy given a vector of prims
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    printf("Inside BVH::build\n");

    // NOTE (PathTracer):
    // This BVH is parameterized on the type of the primitive it contains. This allows
    // us to build a BVH over any type that defines a certain interface. Specifically,
    // we use this to both build a BVH over triangles within each Tri_Mesh, and over
    // a variety of Objects (which might be Tri_Meshes, Spheres, etc.) in Pathtracer.
    //
    // The Primitive interface must implement these two functions:
    //      BBox bbox() const;
    //      Trace hit(const Ray& ray) const;
    // Hence, you may call bbox() and hit() on any value of type Primitive.

    // Keep these two lines of code in your solution. They clear the list of nodes and
    // initialize member variable 'primitives' as a vector of the scene prims
    nodes.clear();
    primitives = std::move(prims);

    // // TODO (PathTracer): Task 3
    // // Modify the code ahead to construct a BVH from the given vector of primitives and maximum leaf
    // // size configuration.
    // //
    // // Please use the SAH as described in class.  We recomment the binned build from lecture.
    // // In general, here is a rough sketch:
    // //
    // //  For each axis X,Y,Z:
    // //     Try possible splits along axis, evaluate SAH for each
    // //  Take minimum cost across all axes.
    // //  Partition primitives into a left and right child group
    // //  Compute left and right child bboxes
    // //  Make the left and right child nodes.
    // //
    // //
    // // While a BVH is conceptually a tree structure, the BVH class uses a single vector (nodes)
    // // to store all the nodes. Therefore, BVH nodes don't contain pointers to child nodes,
    // // but rather the indices of the
    // // child nodes in this array. Hence, to get the child of a node, you have to
    // // look up the child index in this vector (e.g. nodes[node.l]). Similarly,
    // // to create a new node, don't allocate one yourself - use BVH::new_node, which
    // // returns the index of a newly added node.
    // //
    // // As an example of how to make nodes, the starter code below builds a BVH with a
    // // root node that encloses all the primitives and its two descendants at Level 2.
    // // For now, the split is hardcoded such that the first primitive is put in the left
    // // child of the root, and all the other primitives are in the right child.
    // // There are no further descendants.

    // // edge case
    // if(primitives.empty()) {
    //     return;
    // }

    // // compute bounding box for all primitives
    // BBox bb;
    // for(size_t i = 0; i < primitives.size(); ++i) {
    //     bb.enclose(primitives[i].bbox());
    // }

    // // set up root node (root BVH). Notice that it contains all primitives.
    // size_t root_node_addr = new_node();
    // Node& node = nodes[root_node_addr];
    // node.bbox = bb;
    // node.start = 0;
    // node.size = primitives.size();

    // // Create bounding boxes for children
    // BBox split_leftBox;
    // BBox split_rightBox;

    // // compute bbox for left child
    // Primitive& p = primitives[0];
    // BBox pbb = p.bbox();
    // split_leftBox.enclose(pbb);

    // // compute bbox for right child
    // for(size_t i = 1; i < primitives.size(); ++i) {
    //     Primitive& p = primitives[i];
    //     BBox pbb = p.bbox();
    //     split_rightBox.enclose(pbb);
    // }

    // // Note that by construction in this simple example, the primitives are
    // // contiguous as required. But in the students real code, students are
    // // responsible for reorganizing the primitives in the primitives array so that
    // // after a SAH split is computed, the chidren refer to contiguous ranges of primitives.

    // size_t startl = 0;  // starting prim index of left child
    // size_t rangel = 1;  // number of prims in left child
    // size_t startr = startl + rangel;  // starting prim index of right child
    // size_t ranger = primitives.size() - rangel; // number of prims in right child

    // // create child nodes
    // size_t node_addr_l = new_node();
    // size_t node_addr_r = new_node();
    // nodes[root_node_addr].l = node_addr_l;
    // nodes[root_node_addr].r = node_addr_r;

    // nodes[node_addr_l].bbox = split_leftBox;
    // nodes[node_addr_l].start = startl;
    // nodes[node_addr_l].size = rangel;

    // nodes[node_addr_r].bbox = split_rightBox;
    // nodes[node_addr_r].start = startr;
    // nodes[node_addr_r].size = ranger;


    BBox bb;
    for(size_t i = 0; i < primitives.size(); ++i) {
        bb.enclose(primitives[i].bbox());
    }
    // return 0, should be the index of the root node
    size_t root_node_addr = new_node();
    Node& node = nodes[root_node_addr];
    node.bbox = bb;
    node.start = 0;
    node.size = primitives.size();
    // SAH(root_node_addr, max_leaf_size);
    ispc::SAH(root_node_addr, max_leaf_size);
}

template<typename Primitive>
void BVH<Primitive>::SAH(size_t idx, size_t max_leaf_size) {
    if (nodes[idx].size <= max_leaf_size) {
        return;
    }
    // Create bounding boxes for children
    BBox split_leftBox;
    BBox split_rightBox;
    BBox bbox = nodes[idx].bbox;
    int bucket_num = 8;

    int min_cost_axis = 0;
    float min_cost = 0x7fffffff;
    float min_cost_split = 0;


    size_t rangel = 0;
    size_t ranger = 0;
    // can't use buckets[3][bucket_num], will cause bus fault
    for (int i = 0; i < 3; i++) {
        BBox buckets[bucket_num];
        int prim_count[bucket_num];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < bucket_num; j++) {
                prim_count[j] = 0;
            }
        }
        float my_min = bbox.min[i];
        float my_max = bbox.max[i];
        float interval = (my_max - my_min) / (bucket_num + 0.0);
        for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
            Primitive& p = primitives[j];
            BBox pbb = p.bbox();
            
            int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
            bucket_idx = std::clamp(bucket_idx, 0, bucket_num - 1);
            buckets[bucket_idx].enclose(pbb);
            prim_count[bucket_idx]++;
        }
        
        
        for (int j = 0; j < bucket_num - 1; j++) {
            float left_surface_area = 0;
            float right_surface_area = 0;
            float left_prim_count = 0;
            float right_prim_count = 0;

            BBox left_bbox;
            BBox right_bbox;
    
            for (int k = 0; k <= j; k++) {
                left_bbox.enclose(buckets[k]);
                left_prim_count += prim_count[k];
            }
            for (int k = j + 1; k < bucket_num; k++) {
                right_bbox.enclose(buckets[k]);
                right_prim_count += prim_count[k];
            }

            left_surface_area = left_bbox.surface_area();
            right_surface_area = right_bbox.surface_area();

            float total_cost = left_surface_area * left_prim_count + right_surface_area * right_prim_count;
            if (total_cost < min_cost) {
                min_cost = total_cost;
                min_cost_axis = i;
                min_cost_split = my_min + (j + 1) * interval;
                split_leftBox = left_bbox;
                split_rightBox = right_bbox;
                rangel = left_prim_count;
                ranger = right_prim_count;
            }
        }
    }

    // need to reorganize primitives so that the children are contiguous ranges of primitives
    int first = nodes[idx].start;
    for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) 
    {
        Primitive& p = primitives[j];
        BBox pbb = p.bbox();
        if (pbb.center()[min_cost_axis] < min_cost_split)
        {
            std::swap(primitives[j], primitives[first]);
            ++first;
        }
    }

    size_t startl = nodes[idx].start;
    size_t startr = startl + rangel;

    // create child nodes
    size_t node_addr_l = new_node();
    size_t node_addr_r = new_node();
    nodes[idx].l = node_addr_l;
    nodes[idx].r = node_addr_r;

    nodes[node_addr_l].bbox = split_leftBox;
    nodes[node_addr_l].start = startl;
    nodes[node_addr_l].size = rangel;

    nodes[node_addr_r].bbox = split_rightBox;
    nodes[node_addr_r].start = startr;
    nodes[node_addr_r].size = ranger;

    SAH(node_addr_l, max_leaf_size);
    SAH(node_addr_r, max_leaf_size);
}

template<typename Primitive>
Trace BVH<Primitive>::hit(const Ray& ray) const {

    // TODO (PathTracer): Task 3
    // Implement ray - BVH intersection test. A ray intersects
    // with a BVH aggregate if and only if it intersects a primitive in
    // the BVH that is not an aggregate.

    // The starter code simply iterates through all the primitives.
    // Again, remember you can use hit() on any Primitive value.

    // Trace ret;
    // for(const Primitive& prim : primitives) {
    //     Trace hit = prim.hit(ray);
    //     ret = Trace::min(ret, hit);
    // }
    // return ret;
    return find_hit(ray, 0);
}

template<typename Primitive>
Trace BVH<Primitive>::find_hit(const Ray& ray, size_t idx) const {
    Trace ret;
    if (nodes[idx].is_leaf()) {
        for (unsigned long i = nodes[idx].start; i < nodes[idx].start + nodes[idx].size; i++) {
            Trace hit = primitives[i].hit(ray);
            ret = Trace::min(ret, hit);
        }
        return ret;
    } else {
        Vec2 times_l;
        Vec2 times_r;
        // pass reference
        bool hit_l = nodes[nodes[idx].l].bbox.hit(ray, times_l);
        bool hit_r = nodes[nodes[idx].r].bbox.hit(ray, times_r);

        if (hit_l && hit_r) {
            if (times_l.x < times_r.x) {
                Trace ret_l = find_hit(ray, nodes[idx].l);
                Trace ret_r = find_hit(ray, nodes[idx].r);
                return Trace::min(ret_l, ret_r);
            } else {
                Trace ret_r = find_hit(ray, nodes[idx].r);
                Trace ret_l = find_hit(ray, nodes[idx].l);
                return Trace::min(ret_l, ret_r);
            }
        } else if (hit_l) {
            return find_hit(ray, nodes[idx].l);
        } else if (hit_r) {
            return find_hit(ray, nodes[idx].r);
        } else {
            return ret;
        }
    }
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive>
BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive>
bool BVH<Primitive>::Node::is_leaf() const {
    return l == r;
}

template<typename Primitive>
size_t BVH<Primitive>::new_node(BBox box, size_t start, size_t size, size_t l, size_t r) {
    Node n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.l = l;
    n.r = r;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive>
BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive>
std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive>
void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

template<typename Primitive>
size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                                 const Mat4& trans) const {

    std::stack<std::pair<size_t, size_t>> tstack;
    tstack.push({root_idx, 0});
    size_t max_level = 0;

    if(nodes.empty()) return max_level;

    while(!tstack.empty()) {

        auto [idx, lvl] = tstack.top();
        max_level = std::max(max_level, lvl);
        const Node& node = nodes[idx];
        tstack.pop();

        Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
        GL::Lines& add = lvl == level ? active : lines;

        BBox box = node.bbox;
        box.transform(trans);
        Vec3 min = box.min, max = box.max;

        auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

        edge(min, Vec3{max.x, min.y, min.z});
        edge(min, Vec3{min.x, max.y, min.z});
        edge(min, Vec3{min.x, min.y, max.z});
        edge(max, Vec3{min.x, max.y, max.z});
        edge(max, Vec3{max.x, min.y, max.z});
        edge(max, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
        edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
        edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

        if(node.l && node.r) {
            tstack.push({node.l, lvl + 1});
            tstack.push({node.r, lvl + 1});
        } else {
            for(size_t i = node.start; i < node.start + node.size; i++) {
                size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
                max_level = std::max(c, max_level);
            }
        }
    }
    return max_level;
}

} // namespace PT
