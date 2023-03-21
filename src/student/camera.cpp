
#include "../util/camera.h"
#include "../rays/samplers.h"
#include "../util/rand.h"
#include "debug.h"
#include "math.h"

Ray Camera::generate_ray(Vec2 screen_coord) const {

    // TODO (PathTracer): Task 1
    //
    // The input screen_coord is a normalized screen coordinate [0,1]^2
    //
    // You need to transform this 2D point into a 3D position on the sensor plane, which is
    // located one unit away from the pinhole in camera space (aka view space).
    //
    // You'll need to compute this position based on the vertial field of view
    // (vert_fov) of the camera, and the aspect ratio of the output image (aspect_ratio).
    //
    // Tip: compute the ray direction in view space and use
    // the camera space to world space transform (iview) to transform the ray back into world space.

    // tan need to be in radian
    float height = 2.0f * std::tan(vert_fov * 0.5f * M_PI / 180.0f);
    float width = height * aspect_ratio;

    float x = screen_coord[0] * width - width * 0.5f;
    float y = screen_coord[1] * height - height * 0.5f;

    Vec3 dir = Vec3(x, y, -1.0f).unit();
    Vec3 point = Vec3(0, 0, 0);


    // to support multi-spectral rendering, now each ray has its own wavelength lambda
    // when camera generate a ray, the wavelength of the ray could be any float number
    // between 400nm to 690nm. 
    float r = static_cast <float> (rand()) / static_cast <float> (RAND_MAX);
    float lambda = 400.0f + r*(690.0f - 400.0f);
    Ray ray(point, dir,lambda);
    ray.transform(iview);

    return ray;
    // return Ray();
}
