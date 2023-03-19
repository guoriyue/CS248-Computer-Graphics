#define SIMD_WIDTH 8
// can't use class in ISPC
// will get Error: syntax error, unexpected identifier.

struct Vec2 {
    float x, y;
};

struct Vec3 {
    float x, y, z;
};

struct BBox {
    Vec3 min;
    Vec3 max;
};

struct Node {
    BBox bbox;
    size_t start, size, l, r;
};

struct Ray {
    Vec3 point, dir;
};

struct Trace {
    bool hit;
    float distance;
    Vec3 position, normal, origin;
    int material;
};

inline void swap(float& a, float& b) {
    float tmp = a;
    a = b;
    b = tmp;
}

task void bbox_hit_task(const uniform Ray& ray, uniform BBox ispc_bboxs[], uniform Vec2 ispc_times[],
uniform bool ispc_hits[]) {
    Vec2 times = ispc_times[taskIndex];
    BBox bbox = ispc_bboxs[taskIndex];

    float tmin = (bbox.min.x - ray.point.x) / ray.dir.x; 
    float tmax = (bbox.max.x - ray.point.x) / ray.dir.x; 
 
    if (tmin > tmax) {
        swap(tmin, tmax);
    }
 
    float tymin = (bbox.min.y - ray.point.y) / ray.dir.y; 
    float tymax = (bbox.max.y - ray.point.y) / ray.dir.y; 
 
    if (tymin > tymax) {
        swap(tymin, tymax);
    }
 
    if ((tmin > tymax) || (tymin > tmax)) {
        // return false;
        ispc_hits[taskIndex] = false; 
    }

    tmin = max(tmin, tymin);
    tmax = min(tmax, tymax);
 
    float tzmin = (bbox.min.z - ray.point.z) / ray.dir.z; 
    float tzmax = (bbox.max.z - ray.point.z) / ray.dir.z; 
 
    if (tzmin > tzmax) {
        swap(tzmin, tzmax); 
    }
 
    if ((tmin > tzmax) || (tzmin > tmax)) {
        // return false; 
        ispc_hits[taskIndex] = false;
    }
    
    tzmin = max(tzmin, tmin);
    tzmax = min(tzmax, tmax);
    
    times.x = tmin;
    times.y = tmax;
    // return true;
    ispc_hits[taskIndex] = true;
}

// Trace Triangle::hit(const Ray& ray) const {
//     // printf("Triangle::hit\n");
//     // Vertices of triangle - has postion and surface normal
//     // See rays/tri_mesh.h for a description of this struct
    
//     Tri_Mesh_Vert v_0 = vertex_list[v0];
//     Tri_Mesh_Vert v_1 = vertex_list[v1];
//     Tri_Mesh_Vert v_2 = vertex_list[v2];

//     // here just to avoid unused variable warnings, students should remove the following three lines.
//     (void)v_0;
//     (void)v_1;
//     (void)v_2;
    
//     // TODO (PathTracer): Task 2
//     // Intersect this ray with a triangle defined by the above three points.
//     // Intersection should yield a ray t-value, and a hit point (u,v) on the surface of the triangle

//     // You'll need to fill in a "Trace" struct describing information about the hit (or lack of hit)

//     Trace ret;
//     ret.origin = ray.point;
//     ret.hit = false;       // was there an intersection?
//     ret.distance = 0.0f;   // at what distance did the intersection occur?
//     ret.position = Vec3{}; // where was the intersection?
//     ret.normal = Vec3{};   // what was the surface normal at the intersection?
//                            // (this should be interpolated between the three vertex normals)

//     Vec3 e_1 = v_1.position - v_0.position;
//     Vec3 e_2 = v_2.position - v_0.position;
//     Vec3 s = ray.point - v_0.position;
//     auto determinant = [](auto e_1, auto e_2, auto dir) {
//         return dot(cross(e_1, dir), e_2);
//     };

//     // float denominator = dot(cross(e_1, ray.dir), e_2);
//     float denominator = determinant(e_1, e_2, ray.dir);
//     if (std::abs(denominator) < 1e-6) {
//         // parallel to the plane
//         return ret;
//     }
//     Vec3 uvt = Vec3(-1 * determinant(s, ray.dir, e_2), determinant(e_1, s, ray.dir), -1 * determinant(s, e_1, e_2)) / denominator;
    
//     float u = uvt[0];
//     float v = uvt[1];
//     float t = uvt[2];
    

//     if (0 <= u && u <= 1 && 0 <= v && v <= 1 && u + v <= 1 && ray.dist_bounds[0] <= t && t <= ray.dist_bounds[1]) {
//         // printf("u: %f, v: %f, t: %f\n", u, v, t);
//         // printf("ray.dist_bounds[0] %f, ray.dist_bounds[1] %f\n", ray.dist_bounds[0], ray.dist_bounds[1]);
//         ret.hit = true;
//         ret.distance = t;
//         ret.position = ray.point + ray.dir * t;
//         ret.normal = (1 - u - v) * v_0.normal + u * v_1.normal + v * v_2.normal;
//     }

//     return ret;
// }


export void bbox_hit(const uniform Ray& ray, uniform BBox ispc_bboxs[], uniform Vec2 ispc_times[],
uniform bool ispc_hits[]) {
    // launch[0] bbox_hit_task(ray, nodes, times, 0);
    launch[SIMD_WIDTH] bbox_hit_task(ray, ispc_bboxs, ispc_times, ispc_hits);
}




// task void find_hit_task(uniform size_t idx, uniform size_t max_leaf_size, uniform Node nodes[], uniform size_t depth) {
//     if (depth > 3) {
//         return;
//     }
//     launch[0] find_hit_task(idx, max_leaf_size, nodes, depth + 1);
//     launch[1] find_hit_task(idx, max_leaf_size, nodes, depth + 1);
// }

// export void find_hit(uniform size_t idx, uniform size_t max_leaf_size, uniform Node nodes[]) {
//     launch[0] find_hit_task(idx, max_leaf_size, nodes, 0);
// }








// task void SAH_task(uniform size_t idx, uniform size_t max_leaf_size, uniform size_t depth) {
//     if (depth > 10) {
//         return;
//     }
// //     if (nodes[idx].size <= max_leaf_size) {
// //         return;
// //     }
// //     // Create bounding boxes for children
// //     BBox split_leftBox;
// //     BBox split_rightBox;
// //     BBox bbox = nodes[idx].bbox;
// //     int bucket_num = 8;

// //     int min_cost_axis = 0;
// //     float min_cost = 0x7fffffff;
// //     float min_cost_split = 0;


// //     size_t rangel = 0;
// //     size_t ranger = 0;
// //     // can't use buckets[3][bucket_num], will cause bus fault
// //     // foreach (i = 0 ... buf_size3_full) {
// //     for (int i = 0; i < 3; i++) {
// //         BBox buckets[bucket_num];
// //         int prim_count[bucket_num];
// //         for (int i = 0; i < 3; i++) {
// //             for (int j = 0; j < bucket_num; j++) {
// //                 prim_count[j] = 0;
// //             }
// //         }
// //         float my_min = bbox.min[i];
// //         float my_max = bbox.max[i];
// //         float interval = (my_max - my_min) / (bucket_num + 0.0);
// //         for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) {
// //             Primitive& p = primitives[j];
// //             BBox pbb = p.bbox();
            
// //             int bucket_idx = floor((pbb.center()[i] - my_min) / interval);
// //             bucket_idx = clamp(bucket_idx, 0, bucket_num - 1);
// //             buckets[bucket_idx].enclose(pbb);
// //             prim_count[bucket_idx]++;
// //         }
        
        
// //         for (int j = 0; j < bucket_num - 1; j++) {
// //             float left_surface_area = 0;
// //             float right_surface_area = 0;
// //             float left_prim_count = 0;
// //             float right_prim_count = 0;

// //             BBox left_bbox;
// //             BBox right_bbox;
    
// //             for (int k = 0; k <= j; k++) {
// //                 left_bbox.enclose(buckets[k]);
// //                 left_prim_count += prim_count[k];
// //             }
// //             for (int k = j + 1; k < bucket_num; k++) {
// //                 right_bbox.enclose(buckets[k]);
// //                 right_prim_count += prim_count[k];
// //             }

// //             left_surface_area = left_bbox.surface_area();
// //             right_surface_area = right_bbox.surface_area();

// //             float total_cost = left_surface_area * left_prim_count + right_surface_area * right_prim_count;
// //             if (total_cost < min_cost) {
// //                 min_cost = total_cost;
// //                 min_cost_axis = i;
// //                 min_cost_split = my_min + (j + 1) * interval;
// //                 split_leftBox = left_bbox;
// //                 split_rightBox = right_bbox;
// //                 rangel = left_prim_count;
// //                 ranger = right_prim_count;
// //             }
// //         }
// //     }

// //     // need to reorganize primitives so that the children are contiguous ranges of primitives
// //     int first = nodes[idx].start;
// //     for (unsigned long j = nodes[idx].start; j < nodes[idx].start + nodes[idx].size; j++) 
// //     {
// //         Primitive& p = primitives[j];
// //         BBox pbb = p.bbox();
// //         if (pbb.center()[min_cost_axis] < min_cost_split)
// //         {
// //             swap(primitives[j], primitives[first]);
// //             ++first;
// //         }
// //     }

// //     size_t startl = nodes[idx].start;
// //     size_t startr = startl + rangel;

// //     // create child nodes
// //     size_t node_addr_l = new_node();
// //     size_t node_addr_r = new_node();
// //     nodes[idx].l = node_addr_l;
// //     nodes[idx].r = node_addr_r;

// //     nodes[node_addr_l].bbox = split_leftBox;
// //     nodes[node_addr_l].start = startl;
// //     nodes[node_addr_l].size = rangel;

// //     nodes[node_addr_r].bbox = split_rightBox;
// //     nodes[node_addr_r].start = startr;
// //     nodes[node_addr_r].size = ranger;

// //     SAH(node_addr_l, max_leaf_size);
// //     SAH(node_addr_r, max_leaf_size);
//     launch[0] SAH_task(idx, max_leaf_size, depth + 1);
//     launch[1] SAH_task(idx, max_leaf_size, depth + 1);
// }


// // .data()
// // export void build(vector<Primitive>&& prims, size_t max_leaf_size) {
// export void SAH(uniform size_t idx, uniform size_t max_leaf_size) {
//     // nodes.clear();
//     // primitives = move(prims);

//     // BBox bb;
//     // for(size_t i = 0; i < primitives.size(); ++i) {
//     //     bb.enclose(primitives[i].bbox());
//     // }
//     // // return 0, should be the index of the root node
//     // size_t root_node_addr = new_node();
//     // Node& node = nodes[root_node_addr];
//     // node.bbox = bb;
//     // node.start = 0;
//     // node.size = primitives.size();
//     // SAH(root_node_addr, max_leaf_size);
    
    
//     launch[0] SAH_task(idx, max_leaf_size, 0);
// }

