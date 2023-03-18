
#include "../lib/mathlib.h"
#include "debug.h"

bool BBox::hit(const Ray& ray, Vec2& times) const {
    // printf("BBox::hit()\n");

    // TODO (PathTracer):
    // Implement ray - bounding box intersection test
    // If the ray intersected the bounding box within the range given by
    // [times.x,times.y], update times with the new intersection times.

    // return false;

    float tmin = (min.x - ray.point.x) / ray.dir.x; 
    float tmax = (max.x - ray.point.x) / ray.dir.x; 
 
    if (tmin > tmax) {
        std::swap(tmin, tmax);
    }
 
    float tymin = (min.y - ray.point.y) / ray.dir.y; 
    float tymax = (max.y - ray.point.y) / ray.dir.y; 
 
    if (tymin > tymax) {
        std::swap(tymin, tymax);
    }
 
    if ((tmin > tymax) || (tymin > tmax)) {
        return false; 
    }

    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);
 
    float tzmin = (min.z - ray.point.z) / ray.dir.z; 
    float tzmax = (max.z - ray.point.z) / ray.dir.z; 
 
    if (tzmin > tzmax) {
        std::swap(tzmin, tzmax); 
    }
 
    if ((tmin > tzmax) || (tzmin > tmax)) {
        return false; 
    }
    
    tzmin = std::max(tzmin, tmin);
    tzmax = std::min(tzmax, tmax);
    
    times.x = tmin;
    times.y = tmax;
    return true;
}
