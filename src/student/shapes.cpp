
#include "../rays/shapes.h"
#include "debug.h"

namespace PT {

const char* Shape_Type_Names[(int)Shape_Type::count] = {"None", "Sphere"};

BBox Sphere::bbox() const {

    BBox box;
    box.enclose(Vec3(-radius));
    box.enclose(Vec3(radius));
    return box;
}

// Trace Sphere::hit(const Ray& ray) const {

//     // TODO (PathTracer): Task 2
//     // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

//     // If the ray intersects the sphere twice, ret should
//     // represent the first intersection, but remember to respect
//     // ray.dist_bounds! For example, if there are two intersections,
//     // but only the _later_ one is within ray.dist_bounds, you should
//     // return that one!

//     Trace ret;
//     ret.origin = ray.point;
//     ret.hit = false;       // was there an intersection?
//     ret.distance = 0.0f;   // at what distance did the intersection occur?
//     ret.position = Vec3{}; // where was the intersection?
//     ret.normal = Vec3{};   // what was the surface normal at the intersection?

//     float t0;
//     float t1;

//     Vec3 L = bbox().center() - ray.point;
//     float a = dot(ray.dir, ray.dir);
//     float b = 2 * dot(L, ray.dir);
//     float c = dot(L, L) - 1;

//     float discriminant = b * b - 4 * a * c;
//     if (discriminant < 0) {
//         return ret;
//     }
//     else {
//         t0 = (-b - sqrt(discriminant)) / (2 * a);
//         t1 = (-b + sqrt(discriminant)) / (2 * a);
//     }

//     if (ray.dist_bounds[0] <= t0 && t0 <= ray.dist_bounds[1]) {
//         ret.hit = true;
//         ret.distance = t0;
//         ret.position = ray.at(t0);
//         ret.normal = ret.position.unit();
//         return ret;
//     }
//     if (ray.dist_bounds[0] <= t1 && t1 <= ray.dist_bounds[1]) {
//         ret.hit = true;
//         ret.distance = t1;
//         ret.position = ray.at(t1);
//         ret.normal = ret.position.unit();
//         return ret;
//     }

//     return ret;
// }

Trace Sphere::hit(const Ray& ray) const {
    // printf("Sphere::hit()\n");
    // TODO (PathTracer): Task 2
    // Intersect this ray with a sphere of radius Sphere::radius centered at the origin.

    // If the ray intersects the sphere twice, ret should
    // represent the first intersection, but remember to respect
    // ray.dist_bounds! For example, if there are two intersections,
    // but only the _later_ one is within ray.dist_bounds, you should
    // return that one!

    Trace ret;
    ret.origin = ray.point;
    ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?

    //get o and d from the ray
    Vec3 o = ray.point; 
    Vec3 d = ray.dir;

    //check if the term under the square root is positive
    float term = pow(dot(o, d),2) - o.norm_squared() + pow(radius,2);

    if(term < 0) {
        return ret;
    }

    //get the two solutions
    float t1 = -dot(o, d) + sqrt(term);
    float t2 = -dot(o, d) - sqrt(term);

    // true solution is the closer one
    float t = fmin(t1, t2);
    
    // if t1 not in the bound, t should be t2
    if(t1 < ray.dist_bounds.x || t1 > ray.dist_bounds.y) {
        t = t2;
    }

    // if t2 not in the bound, t should be t1
    else if(t2 < ray.dist_bounds.x || t2 > ray.dist_bounds.y) {
        t = t1;
    }

    // if the other t is not in the bound just return
    if (t < ray.dist_bounds.x || t > ray.dist_bounds.y){
        return ret;
    }

    // calculate the hit point
    Vec3 hit_point = o + t * d; 

    ret.hit = true;
    ret.distance = t;
    ret.position = hit_point;
    ret.normal = hit_point.unit();

    return ret;
}

} // namespace PT
