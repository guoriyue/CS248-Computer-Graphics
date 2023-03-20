
// // Parallel BVH Tree
// #define SIMD_WIDTH 4
// // can't use class in ISPC
// // will get Error: syntax error, unexpected identifier.

// struct Vec2 {
//     float x, y;
// };

// struct Vec3 {
//     float x, y, z;
// };

// struct BBox {
//     Vec3 min;
//     Vec3 max;
// };

// struct Node {
//     BBox bbox;
//     size_t start, size;
//     int child[SIMD_WIDTH];
// };

// struct Ray {
//     Vec3 point, dir;
//     Vec2 dist_bounds;
// };

// struct Trace {
//     bool hit;
//     float distance;
//     Vec3 position, normal, origin;
//     int material;
// };

// struct Tri_Mesh_Vert {
//     Vec3 position;
//     Vec3 normal;
// };

// struct Triangle {
//     Tri_Mesh_Vert v_0, v_1, v_2;
// };

// inline void swap(float& a, float& b) {
//     float tmp = a;
//     a = b;
//     b = tmp;
// }

// inline bool is_leaf(Node ispc_node) {
//     for(int i = 0; i < SIMD_WIDTH; i++) {
//         if(ispc_node.child[0] == 0) {
//             continue;
//         } else {
//             return false;
//         }
//     }
//     return true;
// }

// inline Vec3 Vec3_new(float x, float y, float z) {
//     Vec3 v;
//     v.x = 0.0f;
//     v.y = 0.0f;
//     v.z = 0.0f;
//     return v;
// }

// inline Vec3 Vec3_sub(Vec3 a, Vec3 b) {
//     return Vec3_new(a.x - b.x, a.y - b.y, a.z - b.z);
// }
// inline float Vec3_dot(Vec3 a, Vec3 b) {
//     return a.x * b.x + a.y * b.y + a.z * b.z;
// }
// inline Vec3 Vec3_cross(Vec3 a, Vec3 b) {
//     return Vec3_new(a.y * b.z - a.z * b.y, a.z * b.x - a.x * b.z, a.x * b.y - a.y * b.x);
// }
// inline Vec3 Vec3_mul_float(Vec3 a, float b) {
//     return Vec3_new(a.x * b, a.y * b, a.z * b);
// }

// inline Vec3 Vec3_div_float(Vec3 a, float b) {
//     return Vec3_new(a.x / b, a.y / b, a.z / b);
// }

// inline Vec3 Vec3_add(Vec3 a, Vec3 b) {
//     return Vec3_new(a.x + b.x, a.y + b.y, a.z + b.z);
// }

// inline Vec3 Vec3_add_3(Vec3 a, Vec3 b, Vec3 c) {
//     return Vec3_new(a.x + b.x + c.x, a.y + b.y + c.y, a.z + b.z + c.z);
// }

// inline float determinant(Vec3 a, Vec3 b, Vec3 c) {
//     return Vec3_dot(Vec3_cross(a, b), c);
// }

// inline Trace triangle_hit(Triangle ispc_triangle, const uniform Ray& ispc_ray) {
//     Tri_Mesh_Vert v_0 = ispc_triangle.v_0;
//     Tri_Mesh_Vert v_1 = ispc_triangle.v_1;
//     Tri_Mesh_Vert v_2 = ispc_triangle.v_2;
//     // here just to avoid unused variable warnings, students should remove the following three
//     // lines.
//     (void)v_0;
//     (void)v_1;
//     (void)v_2;

//     Trace ret;
//     ret.origin = ispc_ray.point;
//     ret.hit = false;                  // was there an intersection?
//     ret.distance = 0.0f;              // at what distance did the intersection occur?
//     ret.position = Vec3_new(0, 0, 0); // where was the intersection?
//     ret.normal =
//         Vec3_new(0, 0, 0); // what was the surface normal at the intersection?
//                            // (this should be interpolated between the three vertex normals)
//     Vec3 e_1 = Vec3_sub(v_1.position, v_0.position);
//     Vec3 e_2 = Vec3_sub(v_2.position, v_0.position);
//     Vec3 s = Vec3_sub(ispc_ray.point, v_0.position);

//     float denominator = determinant(e_1, e_2, ispc_ray.dir);
//     if(abs(denominator) < 1e-6) {
//         // parallel to the plane
//         return ret;
//     }
//     Vec3 uvt =
//         Vec3_div_float(Vec3_new(-1 * determinant(s, ispc_ray.dir, e_2),
//                                 determinant(e_1, s, ispc_ray.dir), -1 * determinant(s, e_1, e_2)),
//                        denominator);

//     float u = uvt.x;
//     float v = uvt.y;
//     float t = uvt.z;

//     if(0 <= u && u <= 1 && 0 <= v && v <= 1 && u + v <= 1 && ispc_ray.dist_bounds.x <= t &&
//        t <= ispc_ray.dist_bounds.y) {

//         ret.hit = true;
//         ret.distance = t;
//         ret.position = Vec3_add(ispc_ray.point, Vec3_mul_float(ispc_ray.dir, t));
//         ret.normal = Vec3_add_3(Vec3_mul_float(v_0.normal, (1 - u - v)),
//                                 Vec3_mul_float(v_1.normal, u), Vec3_mul_float(v_2.normal, v));
//     }

//     return ret;
// }

// bool bbox_hit(const uniform Ray& ray, BBox bbox, Vec2 times) {

//     float tmin = (bbox.min.x - ray.point.x) / ray.dir.x;
//     float tmax = (bbox.max.x - ray.point.x) / ray.dir.x;

//     if(tmin > tmax) {
//         swap(tmin, tmax);
//     }

//     float tymin = (bbox.min.y - ray.point.y) / ray.dir.y;
//     float tymax = (bbox.max.y - ray.point.y) / ray.dir.y;

//     if(tymin > tymax) {
//         swap(tymin, tymax);
//     }

//     if((tmin > tymax) || (tymin > tmax)) {
//         return false;
//         // ispc_hits[taskIndex] = false;
//     }

//     tmin = max(tmin, tymin);
//     tmax = min(tmax, tymax);

//     float tzmin = (bbox.min.z - ray.point.z) / ray.dir.z;
//     float tzmax = (bbox.max.z - ray.point.z) / ray.dir.z;

//     if(tzmin > tzmax) {
//         swap(tzmin, tzmax);
//     }

//     if((tmin > tzmax) || (tzmin > tmax)) {
//         return false;
//     }

//     tzmin = max(tzmin, tmin);
//     tzmax = min(tzmax, tmax);

//     times.x = tmin;
//     times.y = tmax;
//     return true;
// }

// Trace trace_min(Trace a, Trace b) {
//     if(a.hit && b.hit) {
//         return a.distance < b.distance ? a : b;
//     } else if(a.hit) {
//         return a;
//     } else if(b.hit) {
//         return b;
//     } else {
//         return a;
//     }
// }

// // task
// void find_hit_task(size_t idx, const uniform Ray& ispc_ray, Vec2 ispc_times[], BBox ispc_bboxs[],
//                    bool ispc_hits[], Triangle ispc_triangles[], Trace ispc_ret[],
//                    Node ispc_nodes[]) {
//     if(ispc_hits[idx]) {
//         if(is_leaf(ispc_nodes[idx])) {
//             for(int i = ispc_nodes[idx].start; i < ispc_nodes[idx].start + ispc_nodes[idx].size;
//                 i++) {
//                 Trace hit = triangle_hit(ispc_triangles[i], ispc_ray);
//                 ispc_ret[0] = trace_min(ispc_ret[0], hit);
//             }

//         } else {
//             foreach(i = 0 ... 4) {
//                 ispc_hits[i] = bbox_hit(ispc_ray, ispc_bboxs[i], ispc_times[i]);
//                 find_hit_task(ispc_nodes[idx].child[i], ispc_ray, ispc_times, ispc_bboxs, ispc_hits,
//                               ispc_triangles, ispc_ret, ispc_nodes);
//             }
//         }
//     }
// }

// export void find_hit(uniform size_t idx, const uniform Ray& ispc_ray, Vec2 ispc_times[],
//                      BBox ispc_bboxs[], bool ispc_hits[], Triangle ispc_triangles[],
//                      Trace ispc_ret[], Node ispc_nodes[]) {

//     if(is_leaf(ispc_nodes[idx])) {
//         for(int i = ispc_nodes[idx].start; i < ispc_nodes[idx].start + ispc_nodes[idx].size; i++) {
//             Trace hit = triangle_hit(ispc_triangles[i], ispc_ray);
//             ispc_ret[0] = trace_min(ispc_ret[0], hit);
//         }
//     } else {
//         foreach(i = 0 ... 4) {
//             ispc_hits[i] = bbox_hit(ispc_ray, ispc_bboxs[i], ispc_times[i]);
//             find_hit_task(ispc_nodes[idx].child[i], ispc_ray, ispc_times, ispc_bboxs, ispc_hits,
//                           ispc_triangles, ispc_ret, ispc_nodes);
//         }
//     }
// }




// parallel Bounding Box Collision
#define SIMD_WIDTH 4
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

task void bbox_hit_task(uniform const Ray& ray, uniform BBox ispc_bboxs[],
                        uniform Vec2 ispc_times[], uniform bool ispc_hits[]) {
    Vec2 times = ispc_times[taskIndex];
    BBox bbox = ispc_bboxs[taskIndex];

    float tmin = (bbox.min.x - ray.point.x) / ray.dir.x;
    float tmax = (bbox.max.x - ray.point.x) / ray.dir.x;

    if(tmin > tmax) {
        swap(tmin, tmax);
    }

    float tymin = (bbox.min.y - ray.point.y) / ray.dir.y;
    float tymax = (bbox.max.y - ray.point.y) / ray.dir.y;

    if(tymin > tymax) {
        swap(tymin, tymax);
    }

    if((tmin > tymax) || (tymin > tmax)) {
        ispc_hits[taskIndex] = false;
    }

    tmin = max(tmin, tymin);
    tmax = min(tmax, tymax);

    float tzmin = (bbox.min.z - ray.point.z) / ray.dir.z;
    float tzmax = (bbox.max.z - ray.point.z) / ray.dir.z;

    if(tzmin > tzmax) {
        swap(tzmin, tzmax);
    }

    if((tmin > tzmax) || (tzmin > tmax)) {
        ispc_hits[taskIndex] = false;
    }

    tzmin = max(tzmin, tmin);
    tzmax = min(tzmax, tmax);

    times.x = tmin;
    times.y = tmax;
    ispc_hits[taskIndex] = true;
}

export void bbox_hit(uniform const Ray& ray, uniform BBox ispc_bboxs[], uniform Vec2 ispc_times[],
                     uniform bool ispc_hits[]) {
    launch[SIMD_WIDTH] bbox_hit_task(ray, ispc_bboxs, ispc_times, ispc_hits);
}