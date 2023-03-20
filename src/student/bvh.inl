
#include "../rays/bvh.h"

#include "../rays/ispc_bvh.h"
#include <stack>
#define SIMD_WIDTH 4

namespace PT {

// construct BVH hierarchy given a vector of prims
template<typename Primitive>
void BVH<Primitive>::build(std::vector<Primitive>&& prims, size_t max_leaf_size) {

    nodes.clear();
    primitives = std::move(prims);

    BBox bb;
    for(size_t i = 0; i < primitives.size(); ++i) {
        bb.enclose(primitives[i].bbox());
    }
    // return 0, should be the index of the root node
    size_t root_node_addr = new_node4();
    Node4& node = nodes[root_node_addr];
    // size_t root_node_addr = new_node8();
    // Node8& node = nodes[root_node_addr];
    // size_t root_node_addr = new_node();
    // Node& node = nodes[root_node_addr];

    node.bbox = bb;
    node.start = 0;
    node.size = primitives.size();

    // SAH(root_node_addr, max_leaf_size);

    // SAH4(root_node_addr, max_leaf_size);

    kmeans_root = new Kmeans(2, SIMD_WIDTH, 5, std::vector<size_t>(0, primitives.size()));
    kmeans_root->m_iterations = 2;
    kmeans_root->m_K = SIMD_WIDTH;
    kmeans_root->m_P = 5;
    for(size_t i = 0; i < primitives.size(); ++i) {
        kmeans_root->indexOfPrimitives.push_back(i);
    }
    kmeans_root->cluster = new Cluster[kmeans_root->m_K];
    kmeans_root->children = new Kmeans*[kmeans_root->m_K];
    BBox world;
    for(size_t i = 0; i < primitives.size(); ++i) {
        world.enclose(primitives[i].bbox());
    }

    auto length = [&](Vec3 v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); };
    auto getRandCentroidsOnMesh = [&](size_t k, size_t p, std::vector<size_t> indexOfPrimitives) {
        std::vector<BBox> kCentroids;
        size_t idx_primitive;

        size_t prim_size =
            indexOfPrimitives[indexOfPrimitives.size() - 1] - indexOfPrimitives[0] + 1;
        idx_primitive = rand() % prim_size + indexOfPrimitives[0];
        kCentroids.push_back(primitives[idx_primitive].bbox());

        for(unsigned long i = 1; i < k; ++i) {
            std::vector<BBox> tempP;
            for(unsigned long j = 0; j < p; j++) {
                idx_primitive = rand() % prim_size + indexOfPrimitives[0];
                tempP.push_back(primitives[idx_primitive].bbox());
            }

            int index = 0;
            double maxDistance = -1.0f;

            for(unsigned long k = 0; k < p; k++) {
                for(unsigned long q = 0; q < kCentroids.size(); q++) {
                    BBox bb;
                    bb.enclose(tempP[k]);
                    bb.enclose(kCentroids[q]);
                    double distance = length(bb.max - bb.min);
                    if(distance > maxDistance) {
                        maxDistance = distance;
                        index = k;
                    }
                }
            }
            kCentroids.push_back(tempP[index]);
        }
        return kCentroids;
    };

    std::vector<BBox> kCentroids =
        getRandCentroidsOnMesh(kmeans_root->m_K, kmeans_root->m_P, kmeans_root->indexOfPrimitives);

    for(size_t i = 0; i < kmeans_root->m_K; ++i) {
        kmeans_root->cluster[i].representive = kCentroids[i];
    }
    constructKaryTree(kmeans_root, max_leaf_size);
}

template<typename Primitive>
void BVH<Primitive>::allocateKaryTree(Kmeans* kmeans, size_t max_leaf_size) {
    auto length = [&](Vec3 v) { return sqrt(v.x * v.x + v.y * v.y + v.z * v.z); };
    auto calDistance = [&](BBox b1, BBox b2) {
        double min_value = length(b1.min - b2.min);
        double max_value = length(b1.max - b2.max);
        double res = min_value + max_value;
        return res;
    };

    for(size_t iter = 0; iter < kmeans->m_iterations; ++iter) {
        for(size_t i = 0; i < kmeans->m_K; ++i) {

            if(iter != 0) {
                kmeans->cluster[i].updateRepresentive();
            }
            kmeans->cluster[i].reset();
            kmeans->cluster[i].indexOfPrimitives.clear();
        }

        for(size_t idx_primitives = 0; idx_primitives < kmeans->indexOfPrimitives.size();
            ++idx_primitives) {
            size_t index = 0;
            double minDistance = std::numeric_limits<double>::max();
            BBox temp = primitives[kmeans->indexOfPrimitives[idx_primitives]].bbox();

            for(size_t idx_clusters = 0; idx_clusters < kmeans->m_K; ++idx_clusters) {
                double dist = calDistance(temp, kmeans->cluster[idx_clusters].representive);
                if(dist < minDistance) {
                    minDistance = dist;
                    index = idx_clusters;
                }
            }
            // allocate to correct cluster
            kmeans->cluster[index].add(
                idx_primitives, primitives[kmeans->indexOfPrimitives[idx_primitives]].bbox());
        }
    }
}

template<typename Primitive>
void BVH<Primitive>::constructKaryTree(Kmeans* kmeans, size_t max_leaf_size) {
    allocateKaryTree(kmeans, max_leaf_size);
    for(size_t i = 0; i < kmeans->m_K; i++) {
        if(kmeans->cluster[i].indexOfPrimitives.size() < max_leaf_size * kmeans->m_K) {
            kmeans->children[i] = NULL;
            continue;
        }
        if(kmeans->cluster[i].representive.min.x > 0x3fffff) {
            kmeans->children[i] = NULL;
            continue;
        }

        std::vector<size_t> pTemp;
        for(size_t p = 0; p < kmeans->cluster[i].indexOfPrimitives.size(); ++p) {
            pTemp.push_back(kmeans->cluster[i].indexOfPrimitives[p]);
        }
        kmeans->children[i] = new Kmeans(kmeans->m_iterations, kmeans->m_K, kmeans->m_P, pTemp);

        constructKaryTree(kmeans->children[i], max_leaf_size);
    }
}

template<typename Primitive>
void BVH<Primitive>::SAH(const size_t& idx, const size_t& max_leaf_size) {
    if(nodes[idx].size <= max_leaf_size) {
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
    for(int i = 0; i < 3; i++) {
        BBox buckets[bucket_num];
        int prim_count[bucket_num];
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < bucket_num; j++) {
                prim_count[j] = 0;
            }
        }
        float my_min = bbox.min[i];
        float my_max = bbox.max[i];
        float interval = (my_max - my_min) / (bucket_num + 0.0);
        for(unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
            Primitive& p = primitives[j];
            BBox pbb = p.bbox();

            int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
            bucket_idx = std::clamp(bucket_idx, 0, bucket_num - 1);

            buckets[bucket_idx].enclose(pbb);
            prim_count[bucket_idx]++;
        }

        for(int j = 0; j < bucket_num - 1; j++) {
            float left_surface_area = 0;
            float right_surface_area = 0;
            float left_prim_count = 0;
            float right_prim_count = 0;

            BBox left_bbox;
            BBox right_bbox;

            for(int k = 0; k <= j; k++) {
                left_bbox.enclose(buckets[k]);
                left_prim_count += prim_count[k];
            }
            for(int k = j + 1; k < bucket_num; k++) {
                right_bbox.enclose(buckets[k]);
                right_prim_count += prim_count[k];
            }

            left_surface_area = left_bbox.surface_area();
            right_surface_area = right_bbox.surface_area();

            float total_cost =
                left_surface_area * left_prim_count + right_surface_area * right_prim_count;
            if(total_cost < min_cost) {
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
    for(unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
        Primitive& p = primitives[j];
        BBox pbb = p.bbox();
        if(pbb.center()[min_cost_axis] < min_cost_split) {
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

template<typename Primitive> Trace BVH<Primitive>::hit(const Ray& ray) const {
    // // recursive version
    // return find_hit(ray, 0);

    // // stack version, sequential
    // Trace ret;
    // std::stack<size_t> node_stack;
    // node_stack.push(root_idx);

    // while(!node_stack.empty()) {
    //     auto node_idx = node_stack.top();
    //     node_stack.pop();
    //     auto& node = nodes[node_idx];

    //     Vec2 times0{};
    //     auto hit0 = node.bbox.hit(ray, times0);
    //     if(!hit0 || (ret.hit && ret.distance <= times0.x)) continue;

    //     if(node.is_leaf()) {
    //         auto node_end = node.start + node.size;
    //         for(auto i = node.start; i < node_end; ++i) {
    //             Trace hit = primitives[i].hit(ray);
    //             ret = Trace::min(ret, hit);
    //         }
    //     } else {
    //         Vec2 times1, times2;
    //         auto hit1 = nodes[node.l].bbox.hit(ray, times1);
    //         auto hit2 = nodes[node.r].bbox.hit(ray, times2);

    //         if(hit1 && hit2) {
    //             auto first = times1.x < times2.x ? node.l : node.r;
    //             auto second = times1.x < times2.x ? node.r : node.l;
    //             node_stack.push(first);
    //             node_stack.push(second);
    //         } else if(hit1) {
    //             node_stack.push(node.l);
    //         } else if(hit2) {
    //             node_stack.push(node.r);
    //         }
    //     }
    // }
    // return ret;

    return hit_kmeans(ray, kmeans_root);
    // return hit_queue4(ray);
    // return hit_queue8(ray);
}

template<typename Primitive>
Trace BVH<Primitive>::hit_kmeans(const Ray& ray, Kmeans* kmeans) const {
    // // sequential version
    // Trace ret;
    // bool visit = false;
    // for(int j = 0; j < SIMD_WIDTH; j++) {
    //     if(kmeans->children[j] != NULL && kmeans->cluster[j].indexOfPrimitives.size()) {
    //         visit = true;

    //         Vec2 times;
    //         bool hit = kmeans->cluster[j].representive.hit(ray, times);
    //         if(hit) {
    //             Trace hit = hit_kmeans(ray, kmeans->children[j]);
    //             ret = Trace::min(ret, hit);
    //         }
    //     }
    // }
    // if(!visit)
    // {
    //     for(size_t p=0;p<kmeans->indexOfPrimitives.size();++p){
    //         Trace hit = primitives[kmeans->indexOfPrimitives[p]].hit(ray);
    //         ret = Trace::min(ret, hit);
    //     }
    // }
    // return ret;

    auto ispc_Vec3 = [](const Vec3& v) {
        ispc::Vec3 res;
        res.x = v.x;
        res.y = v.y;
        res.z = v.z;
        return res;
    };

    auto is_leaf = [](Kmeans* kmeans) {
        return kmeans->children[0] == NULL && kmeans->children[1] == NULL &&
               kmeans->children[2] == NULL && kmeans->children[3] == NULL;
    };

    ispc::Ray ispc_ray;
    ispc::Vec2* ispc_times = new ispc::Vec2[SIMD_WIDTH];
    bool* ispc_hits = new bool[SIMD_WIDTH];
    Vec2* times = new Vec2[SIMD_WIDTH];

    ispc_ray.point = ispc_Vec3(ray.point);
    ispc_ray.dir = ispc_Vec3(ray.dir);
    Trace ret;
    std::stack<Kmeans*> node_stack;
    node_stack.push(kmeans);
    while(!node_stack.empty()) {
        Kmeans* Kmeans_top = node_stack.top();
        node_stack.pop();
        if(is_leaf(Kmeans_top)) {

            for(size_t p = 0; p < Kmeans_top->indexOfPrimitives.size(); ++p) {
                Trace hit = primitives[Kmeans_top->indexOfPrimitives[p]].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            ispc::BBox* ispc_bbox = new ispc::BBox[SIMD_WIDTH];
            std::vector<ispc::BBox> ispc_bbox4;
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(kmeans->children[i] && kmeans->cluster[i].indexOfPrimitives.size()) {
                    ispc_bbox[i].min = ispc_Vec3(kmeans->cluster[i].representive.min);
                    ispc_bbox[i].max = ispc_Vec3(kmeans->cluster[i].representive.max);
                }
            }
            ispc::bbox_hit(ispc_ray, ispc_bbox, ispc_times, ispc_hits);

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                times[i].x = ispc_times[i].x;
                times[i].y = ispc_times[i].y;
            }

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(ispc_hits[i]) {
                    if(kmeans->children[i]) {
                        node_stack.push(kmeans->children[i]);
                    }
                }
            }
        }
    }
    return ret;
}

template<typename Primitive>
Trace BVH<Primitive>::find_hit(const Ray& ray, const size_t& idx) const {
    Trace ret;
    if(nodes[idx].is_leaf()) {
        for(unsigned long i = nodes[idx].start; i < nodes[idx].start + nodes[idx].size; i++) {
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

        if(hit_l && hit_r) {
            if(times_l.x < times_r.x) {
                Trace ret_l = find_hit(ray, nodes[idx].l);
                Trace ret_r = find_hit(ray, nodes[idx].r);
                return Trace::min(ret_l, ret_r);
            } else {
                Trace ret_r = find_hit(ray, nodes[idx].r);
                Trace ret_l = find_hit(ray, nodes[idx].l);
                return Trace::min(ret_l, ret_r);
            }
        } else if(hit_l) {
            return find_hit(ray, nodes[idx].l);
        } else if(hit_r) {
            return find_hit(ray, nodes[idx].r);
        } else {
            return ret;
        }
    }
}

template<typename Primitive> Trace BVH<Primitive>::hit_queue(const Ray& ray) const {
    Trace ret;
    std::stack<size_t> node_stack;
    node_stack.push(0);
    while(!node_stack.empty()) {
        size_t node_idx = node_stack.top();
        node_stack.pop();
        const Node& node = nodes[node_idx];

        // with early return
        Vec2 times0{};
        bool hit0 = node.bbox.hit(ray, times0);
        if(!hit0 || (ret.hit && ret.distance <= times0.x)) {
            continue;
        }

        if(node.is_leaf()) {
            size_t node_end = node.start + node.size;
            for(size_t i = node.start; i < node_end; ++i) {
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            Vec2 times1, times2;
            bool hit1 = nodes[node.l].bbox.hit(ray, times1);
            bool hit2 = nodes[node.r].bbox.hit(ray, times2);
            if(hit1 && hit2) {
                size_t first = times1.x < times2.x ? node.l : node.r;
                size_t second = times1.x < times2.x ? node.r : node.l;
                node_stack.push(first);
                node_stack.push(second);
            } else if(hit1) {
                node_stack.push(node.l);
            } else if(hit2) {
                node_stack.push(node.r);
            }
        }
    }
    return ret;
}

template<typename Primitive>
BVH<Primitive>::BVH(std::vector<Primitive>&& prims, size_t max_leaf_size) {
    build(std::move(prims), max_leaf_size);
}

template<typename Primitive> BVH<Primitive> BVH<Primitive>::copy() const {
    BVH<Primitive> ret;
    ret.nodes = nodes;
    ret.primitives = primitives;
    ret.root_idx = root_idx;
    return ret;
}

template<typename Primitive> bool BVH<Primitive>::Node::is_leaf() const {
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

template<typename Primitive> BBox BVH<Primitive>::bbox() const {
    return nodes[root_idx].bbox;
}

template<typename Primitive> std::vector<Primitive> BVH<Primitive>::destructure() {
    nodes.clear();
    return std::move(primitives);
}

template<typename Primitive> void BVH<Primitive>::clear() {
    nodes.clear();
    primitives.clear();
}

// // original visualize function
// template<typename Primitive>
// size_t BVH<Primitive>::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
//                                  const Mat4& trans) const {

//     std::stack<std::pair<size_t, size_t>> tstack;
//     tstack.push({root_idx, 0});
//     size_t max_level = 0;

//     if(nodes.empty()) return max_level;

//     while(!tstack.empty()) {

//         auto [idx, lvl] = tstack.top();
//         max_level = std::max(max_level, lvl);
//         const Node& node = nodes[idx];
//         tstack.pop();

//         Vec3 color = lvl == level ? Vec3(1.0f, 0.0f, 0.0f) : Vec3(1.0f);
//         GL::Lines& add = lvl == level ? active : lines;

//         BBox box = node.bbox;
//         box.transform(trans);
//         Vec3 min = box.min, max = box.max;

//         auto edge = [&](Vec3 a, Vec3 b) { add.add(a, b, color); };

//         edge(min, Vec3{max.x, min.y, min.z});
//         edge(min, Vec3{min.x, max.y, min.z});
//         edge(min, Vec3{min.x, min.y, max.z});
//         edge(max, Vec3{min.x, max.y, max.z});
//         edge(max, Vec3{max.x, min.y, max.z});
//         edge(max, Vec3{max.x, max.y, min.z});
//         edge(Vec3{min.x, max.y, min.z}, Vec3{max.x, max.y, min.z});
//         edge(Vec3{min.x, max.y, min.z}, Vec3{min.x, max.y, max.z});
//         edge(Vec3{min.x, min.y, max.z}, Vec3{max.x, min.y, max.z});
//         edge(Vec3{min.x, min.y, max.z}, Vec3{min.x, max.y, max.z});
//         edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, max.y, min.z});
//         edge(Vec3{max.x, min.y, min.z}, Vec3{max.x, min.y, max.z});

//         if(node.l && node.r) {
//             tstack.push({node.l, lvl + 1});
//             tstack.push({node.r, lvl + 1});
//         } else {
//             for(size_t i = node.start; i < node.start + node.size; i++) {
//                 size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
//                 max_level = std::max(c, max_level);
//             }
//         }
//     }
//     return max_level;
// }

// visualize function for node4
// does not actually visualize primitives, just to avoid errors
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
        const Node4& node = nodes[idx];
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

        // if((node.child[0] && (!node.child[1] && !node.child[2] && !node.child[3]))
        // || (node.child[1] && (!node.child[0] && !node.child[2] && !node.child[3]))
        // || (node.child[2] && (!node.child[0] && !node.child[1] && !node.child[3]))
        // || (node.child[3] && (!node.child[0] && !node.child[1] && !node.child[2]))) {
        //     for(size_t i = node.start; i < node.start + node.size; i++) {
        //         size_t c = primitives[i].visualize(lines, active, level - lvl, trans);
        //         max_level = std::max(c, max_level);
        //     }
        // }
        // else
        // {
        //     for(size_t i = 0; i < SIMD_WIDTH; i++) {
        //         if(node.child[i]) {
        //             tstack.push({node.child[i], lvl + 1});
        //         }
        //     }
        // }
    }
    return max_level;
}

template<typename Primitive> bool BVH<Primitive>::Node4::is_leaf4() const {
    return child[0] == child[1] && child[1] == child[2] && child[2] == child[3];
}

template<typename Primitive> size_t BVH<Primitive>::new_node4(BBox box, size_t start, size_t size) {
    Node4 n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    n.child[0] = 0;
    n.child[1] = 0;
    n.child[2] = 0;
    n.child[3] = 0;
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive> Trace BVH<Primitive>::hit_queue4(const Ray& ray) const {

    auto ispc_Vec3 = [](const Vec3& v) {
        ispc::Vec3 res;
        res.x = v.x;
        res.y = v.y;
        res.z = v.z;
        return res;
    };

    // // for parallel whold tree
    // auto ispc_Vec2 = [](const Vec2& v) {
    //     ispc::Vec2 res;
    //     res.x = v.x; res.y = v.y;
    //     return res;
    // };

    // ispc::Ray ispc_ray;
    // ispc::Vec2 *ispc_times = new ispc::Vec2[SIMD_WIDTH];
    // ispc::BBox *ispc_bbox = new ispc::BBox[SIMD_WIDTH];
    // bool* ispc_hits = new bool[SIMD_WIDTH];
    // Vec2 *times = new Vec2[SIMD_WIDTH];
    // ispc::Trace *ispc_ret = new ispc::Trace[1];
    // ispc::Triangle* ispc_triangles = new ispc::Triangle[primitives.size()];
    // ispc::Node* ispc_nodes = new ispc::Node[nodes.size()];
    // for(int i = 0; i < primitives.size(); i++) {
    //     ispc_triangles[i].v_0 = vertex_list[primitives[i].v0];
    //     ispc_triangles[i].v_1 = vertex_list[primitives[i].v1];
    //     ispc_triangles[i].v_2 = vertex_list[primitives[i].v2];
    // }
    // for(int i = 0; i < nodes.size(); i++) {
    //     ispc_nodes[i].bbox.min = ispc_Vec3(nodes[i].bbox.min);
    //     ispc_nodes[i].bbox.max = ispc_Vec3(nodes[i].bbox.max);
    //     ispc_nodes[i].start = nodes[i].start;
    //     ispc_nodes[i].size = nodes[i].size;
    //     for(int j = 0; j < SIMD_WIDTH; ++j) {
    //         ispc_nodes[i].child[j] = nodes[i].child[j];
    //     }
    // }
    // ispc_ray.point = ispc_Vec3(ray.point);
    // ispc_ray.dir = ispc_Vec3(ray.dir);
    // ispc_ray.dist_bounds = ispc_Vec2(ray.dist_bounds);

    // for(int i = 0; i < SIMD_WIDTH; ++i) {
    //     ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
    //     ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
    // }

    // ispc::find_hit(0, ispc_ray, ispc_times, ispc_bbox, ispc_hits, ispc_triangles, ispc_ret,
    // ispc_nodes); Trace ret; ret.hit = ispc_ret[0].hit; ret.distance = ispc_ret[0].distance;
    // ret.position = Vec3(ispc_ret[0].position.x, ispc_ret[0].position.y, ispc_ret[0].position.z);
    // ret.normal = Vec3(ispc_ret[0].normal.x, ispc_ret[0].normal.y, ispc_ret[0].normal.z);
    // ret.origin = Vec3(ispc_ret[0].origin.x, ispc_ret[0].origin.y, ispc_ret[0].origin.z);
    // ret.material = ispc_ret[0].material;
    // return ret;

    ispc::Ray ispc_ray;
    ispc::Vec2* ispc_times = new ispc::Vec2[SIMD_WIDTH];
    ispc::BBox* ispc_bbox = new ispc::BBox[nodes.size()];
    for(unsigned long i = 0; i < nodes.size(); i++) {
        ispc_bbox[i].min = ispc_Vec3(nodes[i].bbox.min);
        ispc_bbox[i].max = ispc_Vec3(nodes[i].bbox.max);
    }
    bool* ispc_hits = new bool[SIMD_WIDTH];
    Vec2* times = new Vec2[SIMD_WIDTH];

    ispc_ray.point = ispc_Vec3(ray.point);
    ispc_ray.dir = ispc_Vec3(ray.dir);
    Trace ret;
    std::stack<size_t> node_stack;
    node_stack.push(0);
    while(!node_stack.empty()) {
        size_t node_idx = node_stack.top();
        node_stack.pop();
        const Node4& node = nodes[node_idx];
        // with early return
        Vec2 times0{};
        bool hit0 = node.bbox.hit(ray, times0);
        if(!hit0 || (ret.hit && ret.distance <= times0.x)) {
            continue;
        }
        if(node.is_leaf4()) {
            size_t node_end = node.start + node.size;
            for(size_t i = node.start; i < node_end; ++i) {
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            std::vector<ispc::BBox> ispc_bbox4;
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(node.child[i]) {
                    ispc_bbox4.push_back(ispc_bbox[node.child[i]]);
                }
            }
            ispc::bbox_hit(ispc_ray, ispc_bbox4.data(), ispc_times, ispc_hits);
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                times[i].x = ispc_times[i].x;
                times[i].y = ispc_times[i].y;
            }

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(ispc_hits[i]) {
                    if(node.child[i]) {
                        node_stack.push(node.child[i]);
                    }
                }
            }
        }
    }
    return ret;
}

template<typename Primitive>
void BVH<Primitive>::SAH4(const size_t& idx, const size_t& max_leaf_size) {
    if(nodes[idx].size <= max_leaf_size) {
        return;
    }

    // Create bounding boxes for children
    BBox* split_bbox = new BBox[SIMD_WIDTH];
    size_t* rangel = new size_t[SIMD_WIDTH];
    size_t* ranger = new size_t[SIMD_WIDTH];
    size_t* node_addr = new size_t[SIMD_WIDTH];
    for(int i = 0; i < SIMD_WIDTH; ++i) {
        node_addr[i] = new_node4();
    }

    size_t* startl = new size_t[SIMD_WIDTH];
    size_t* startr = new size_t[SIMD_WIDTH];
    bucket_split(idx, rangel[0], ranger[0], split_bbox[0], split_bbox[1]);

    startl[0] = nodes[idx].start;
    startr[0] = startl[0] + rangel[0];

    for(int i = 0; i < 2; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i % 2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i / 2))];
            nodes[node_addr[i]].size = rangel[int(floor(i / 2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i / 2))];
            nodes[node_addr[i]].size = ranger[int(floor(i / 2))];
        }
    }

    bucket_split(node_addr[0], rangel[0], ranger[0], split_bbox[0], split_bbox[1]);
    bucket_split(node_addr[1], rangel[1], ranger[1], split_bbox[2], split_bbox[3]);
    for(int i = 0; i < 2; ++i) {
        startl[i] = nodes[node_addr[i]].start;
        startr[i] = startl[i] + rangel[i];
    }
    for(int i = 0; i < 4; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i % 2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i / 2))];
            nodes[node_addr[i]].size = rangel[int(floor(i / 2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i / 2))];
            nodes[node_addr[i]].size = ranger[int(floor(i / 2))];
        }
    }

    for(int i = 0; i < 4; ++i) {
        if(nodes[node_addr[i]].size) {
            nodes[idx].child[i] = node_addr[i];
        }
    }

    for(int i = 0; i < 4; ++i) {
        if(nodes[node_addr[i]].size) {
            SAH4(nodes[idx].child[i], max_leaf_size);
        }
    }
}

template<typename Primitive> bool BVH<Primitive>::Node8::is_leaf8() const {
    return child[0] == 0 && child[0] == child[1] && child[1] == child[2] && child[2] == child[3] &&
           child[3] == child[4] && child[4] == child[5] && child[5] == child[6] &&
           child[6] == child[7];
}

template<typename Primitive> size_t BVH<Primitive>::new_node8(BBox box, size_t start, size_t size) {
    Node8 n;
    n.bbox = box;
    n.start = start;
    n.size = size;
    for(size_t i = 0; i < SIMD_WIDTH; i++) {
        n.child[i] = 0;
    }
    nodes.push_back(n);
    return nodes.size() - 1;
}

template<typename Primitive> Trace BVH<Primitive>::hit_queue8(const Ray& ray) const {
    auto ispc_Vec3 = [](const Vec3& v) {
        ispc::Vec3 res;
        res.x = v.x;
        res.y = v.y;
        res.z = v.z;
        return res;
    };

    ispc::Vec2* ispc_times = new ispc::Vec2[SIMD_WIDTH];
    ispc::BBox* ispc_bbox = new ispc::BBox[SIMD_WIDTH];
    bool* ispc_hits = new bool[SIMD_WIDTH];
    Vec2* times = new Vec2[SIMD_WIDTH];
    ispc::Ray ispc_ray;

    Trace ret;
    std::stack<size_t> node_stack;
    node_stack.push(0);
    while(!node_stack.empty()) {
        size_t node_idx = node_stack.top();
        node_stack.pop();
        const Node8& node = nodes[node_idx];
        // with early return
        Vec2 times0{};
        bool hit0 = node.bbox.hit(ray, times0);
        if(!hit0 || (ret.hit && ret.distance <= times0.x)) {
            continue;
        }
        if(node.is_leaf8()) {
            size_t node_end = node.start + node.size;
            for(size_t i = node.start; i < node_end; ++i) {
                Trace hit = primitives[i].hit(ray);
                ret = Trace::min(ret, hit);
            }
        } else {
            ispc_ray.point = ispc_Vec3(ray.point);
            ispc_ray.dir = ispc_Vec3(ray.dir);

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                ispc_bbox[i].min = ispc_Vec3(nodes[node.child[i]].bbox.min);
                ispc_bbox[i].max = ispc_Vec3(nodes[node.child[i]].bbox.max);
            }

            ispc::bbox_hit(ispc_ray, ispc_bbox, ispc_times, ispc_hits);
            for(int i = 0; i < SIMD_WIDTH; ++i) {
                times[i].x = ispc_times[i].x;
                times[i].y = ispc_times[i].y;
            }

            for(int i = 0; i < SIMD_WIDTH; ++i) {
                if(ispc_hits[i]) {
                    if(node.child[i]) {
                        node_stack.push(node.child[i]);
                    }
                }
            }
        }
    }
    return ret;
}

template<typename Primitive>
void BVH<Primitive>::SAH8(const size_t& idx, const size_t& max_leaf_size) {
    if(nodes[idx].size <= max_leaf_size) {
        return;
    }

    // Create bounding boxes for children
    BBox* split_bbox = new BBox[SIMD_WIDTH];
    size_t* rangel = new size_t[SIMD_WIDTH];
    size_t* ranger = new size_t[SIMD_WIDTH];
    size_t* node_addr = new size_t[SIMD_WIDTH];
    for(int i = 0; i < SIMD_WIDTH; ++i) {
        node_addr[i] = new_node8();
    }

    size_t* startl = new size_t[SIMD_WIDTH];
    size_t* startr = new size_t[SIMD_WIDTH];
    bucket_split(idx, rangel[0], ranger[0], split_bbox[0], split_bbox[1]);

    startl[0] = nodes[idx].start;
    startr[0] = startl[0] + rangel[0];

    for(int i = 0; i < 2; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i % 2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i / 2))];
            nodes[node_addr[i]].size = rangel[int(floor(i / 2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i / 2))];
            nodes[node_addr[i]].size = ranger[int(floor(i / 2))];
        }
    }

    bucket_split(node_addr[0], rangel[0], ranger[0], split_bbox[0], split_bbox[1]);
    bucket_split(node_addr[1], rangel[1], ranger[1], split_bbox[2], split_bbox[3]);
    for(int i = 0; i < 2; ++i) {
        startl[i] = nodes[node_addr[i]].start;
        startr[i] = startl[i] + rangel[i];
    }
    for(int i = 0; i < 4; ++i) {
        nodes[node_addr[i]].bbox = split_bbox[i];
        if(i % 2 == 0) {
            nodes[node_addr[i]].start = startl[int(floor(i / 2))];
            nodes[node_addr[i]].size = rangel[int(floor(i / 2))];
        } else {
            nodes[node_addr[i]].start = startr[int(floor(i / 2))];
            nodes[node_addr[i]].size = ranger[int(floor(i / 2))];
        }
    }

    bucket_split(node_addr[0], rangel[0], ranger[0], split_bbox[0], split_bbox[1]);
    bucket_split(node_addr[1], rangel[1], ranger[1], split_bbox[2], split_bbox[3]);
    bucket_split(node_addr[2], rangel[2], ranger[2], split_bbox[4], split_bbox[5]);
    bucket_split(node_addr[3], rangel[3], ranger[3], split_bbox[6], split_bbox[7]);
    for(int i = 0; i < 4; ++i) {
        // printf("here %d\n", i);
        startl[i] = nodes[node_addr[i]].start;
        startr[i] = startl[i] + rangel[i];
    }

    // printf("reach here\n");
    for(int i = 0; i < 8; ++i) {
        // printf("node_addr[i] = %zu\n", node_addr[i]);
        if(node_addr[i] < 0xffffffff) {
            nodes[node_addr[i]].bbox = split_bbox[i];
            if(i % 2 == 0) {
                nodes[node_addr[i]].start = startl[int(floor(i / 2))];
                nodes[node_addr[i]].size = rangel[int(floor(i / 2))];
            } else {
                nodes[node_addr[i]].start = startr[int(floor(i / 2))];
                nodes[node_addr[i]].size = ranger[int(floor(i / 2))];
            }
        }
    }

    for(int i = 0; i < 8; ++i) {
        if(nodes[node_addr[i]].size) {
            nodes[idx].child[i] = node_addr[i];
        }
    }

    for(int i = 0; i < 8; ++i) {
        if(nodes[node_addr[i]].size) {
            SAH8(nodes[idx].child[i], max_leaf_size);
        }
    }
}

template<typename Primitive>
void BVH<Primitive>::bucket_split(const size_t& idx, size_t& rangel, size_t& ranger,
                                  BBox& split_leftBox, BBox& split_rightBox) {
    int bucket_num = 8;
    BBox bbox = nodes[idx].bbox;
    int min_cost_axis = 0;
    float min_cost = 0x7fffffff;
    float min_cost_split = 0;
    // can't use buckets[3][bucket_num], will cause bus fault
    for(int i = 0; i < 3; i++) {
        BBox buckets[bucket_num];
        int prim_count[bucket_num];
        for(int i = 0; i < 3; i++) {
            for(int j = 0; j < bucket_num; j++) {
                prim_count[j] = 0;
            }
        }
        float my_min = bbox.min[i];
        float my_max = bbox.max[i];
        float interval = (my_max - my_min) / (bucket_num + 0.0);
        for(unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
            Primitive& p = primitives[j];
            BBox pbb = p.bbox();

            int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
            bucket_idx = std::clamp(bucket_idx, 0, bucket_num - 1);
            buckets[bucket_idx].enclose(pbb);
            prim_count[bucket_idx]++;
        }

        for(int j = 0; j < bucket_num - 1; j++) {
            float left_surface_area = 0;
            float right_surface_area = 0;
            float left_prim_count = 0;
            float right_prim_count = 0;

            BBox left_bbox;
            BBox right_bbox;

            for(int k = 0; k <= j; k++) {
                left_bbox.enclose(buckets[k]);
                left_prim_count += prim_count[k];
            }
            for(int k = j + 1; k < bucket_num; k++) {
                right_bbox.enclose(buckets[k]);
                right_prim_count += prim_count[k];
            }

            left_surface_area = left_bbox.surface_area();
            right_surface_area = right_bbox.surface_area();

            float total_cost =
                left_surface_area * left_prim_count + right_surface_area * right_prim_count;
            if(total_cost < min_cost) {
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
    for(unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
        Primitive& p = primitives[j];
        BBox pbb = p.bbox();
        if(pbb.center()[min_cost_axis] < min_cost_split) {
            std::swap(primitives[j], primitives[first]);
            ++first;
        }
    }
}

} // namespace PT
