
#include "../rays/pathtracer.h"
#include "../rays/samplers.h"
#include "../util/rand.h"
#include "../lib/newspectrum.h"
#include "debug.h"

namespace PT {

// Return the radiance along a ray entering the camera and landing on a
// point within pixel (x,y) of the output image.
//
RT_Result Pathtracer::trace_pixel(size_t x, size_t y) {

    Vec2 xy((float)x, (float)y);
    Vec2 wh((float)out_w, (float)out_h);

    // TODO (PathTracer): Task 1

    // Generate a sample within the pixel with coordinates xy and return the
    // incoming light using trace_ray.

    // If n_samples is 1, please send the ray through the center of the pixel.
    // If n_samples > 1, please send the ray through any random point within the pixel

    // Tip: consider making a call to Samplers::Rect::Uniform

    // Tip: you may want to use log_ray for debugging. Given ray t, the following lines
    // of code will log .03% of all rays (see util/rand.h) for visualization in the app.
    // see student/debug.h for more detail.
    // if (RNG::coin_flip(0.0003f))
    //    log_ray(out, 10.0f);

    // Spectrum sum_ray;
    // Samplers::Rect::Uniform sampler;
    //n_samples = 1;
    // float pdf = 1.0f;
    // for(size_t i = 0; i < n_samples; i++) {
    //     Vec2 xy((float)x, (float)y);
    //     Vec2 wh((float)out_w, (float)out_h);
    //     if (n_samples == 1) {
    //         xy += Vec2(0.5f, 0.5f);
    //     }
    //     else {
    //         // This is equivalent to saying that the ray tracer wil shoot n_samples camera rays per pixel.
    //         xy += sampler.sample(pdf);
    //     }
    //     Ray out = camera.generate_ray(xy / wh);
    //     // // for shorter time
    //     // if (RNG::coin_flip(0.000003f)) {
    //     //     log_ray(out, 10.0f);
    //     // }
    //     // if (x == 0 && y == 0) {
    //     //     log_ray(out, 10.0f);
    //     // }
    //     sum_ray += trace_ray(out);
    // }
    // return sum_ray / n_samples;

    Samplers::Rect::Uniform sampler;
    float pdf;

    if (n_samples == 1) {
        xy += Vec2(0.5f, 0.5f);
    }
    else {
        xy += sampler.sample(pdf);
    }
    // As an example, the code below generates a ray through the bottom left of the
    // specified pixel
    Ray out = camera.generate_ray(xy / wh);
    // if (RNG::coin_flip(0.0005f)) {
    //     log_ray(out, 10.0f);
    // }
    float p = trace_ray(out);
    //return both ray wavelength and intensity to do_trace function 
    RT_Result ret;
    ret.lambda = out.lambda;
    ret.p = p*2.5;
    return ret;


}
// Since each ray only has one wavelength, trace_ray now returns a float. 
// The number is the light intensity at that particular wavelength.
float Pathtracer::trace_ray(const Ray& ray) {
    // Trace ray into scene. If nothing is hit, sample the environment
    Trace hit = scene.hit(ray);
    if(!hit.hit) {
        if(env_light.has_value()) {
            return Old2NewSpectrum(env_light.value().sample_direction(ray.dir)).sampleAtLambda(ray.lambda);
        }
        return {};
    }

    // If we're using a two-sided material, treat back-faces the same as front-faces
    const BSDF& bsdf = materials[hit.material];
    if(!bsdf.is_sided() && dot(hit.normal, ray.dir) > 0.0f) {
        hit.normal = -hit.normal;
    }

    // Set up a coordinate frame at the hit point, where the surface normal becomes {0, 1, 0}
    // This gives us out_dir and later in_dir in object space, where computations involving the
    // normal become much easier. For example, cos(theta) = dot(N,dir) = dir.y!
    Mat4 object_to_world = Mat4::rotate_to(hit.normal);
    Mat4 world_to_object = object_to_world.T();
    Vec3 out_dir = world_to_object.rotate(ray.point - hit.position).unit();

    // Debugging: if the normal colors flag is set, return the normal color
    if(debug_data.normal_colors) return Old2NewSpectrum(Spectrum::direction(hit.normal)).sampleAtLambda(ray.lambda);

    // Now we can compute the rendering equation at this point.
    // We split it into two stages:
    //  1. sampling direct lighting (i.e. directly connecting the current path to
    //     each light in the scene)
    //  2. sampling the BSDF to create a new path segment

    // TODO (PathTracer): Task 4
    // The starter code sets radiance_out to (0.25,0.25,0.25) so that you can test your geometry
    // queries before you implement real lighting in Tasks 4 and 5. (i.e, anything that gets hit is not black.)
    // You should change this to (0,0,0) and accumulate the direct and indirect lighting computed below.
    float radiance_out = 0.0f;
    {

        // lambda function to sample a light. Called in loop below.
        auto sample_light = [&](const auto& light) {
            // If the light is discrete (e.g. a point light), then we only need
            // one sample, as all samples will be equivalent
            int samples = light.is_discrete() ? 1 : (int)n_area_samples;
            for(int i = 0; i < samples; i++) {

                // Grab a sample of the light source. See rays/light.h for definition of this struct.
                // Most importantly for Task 4, it contains the distance to the light from hit.position. 
                Light_Sample sample = light.sample(hit.position);
                Vec3 in_dir = world_to_object.rotate(sample.direction);

                // If the light is below the horizon, ignore it
                float cos_theta = in_dir.y;
                if(cos_theta <= 0.0f) continue;

                // If the BSDF has 0 throughput in this direction, ignore it.
                // This is another opportunity to do Russian roulette on low-throughput rays,
                // which would allow us to skip the shadow ray cast, increasing efficiency.
                Spectrum attenuation = bsdf.evaluate(out_dir, in_dir);
                if(attenuation.luma() == 0.0f) continue;

                // TODO (PathTracer): Task 4
                // Construct a shadow ray and compute whether the intersected surface is
                // in shadow. Only accumulate light if not in shadow.

                // Tip: since you're creating the shadow ray at the intersection point, it may
                // intersect the surface at time=0. Similarly, if the ray is allowed to have
                // arbitrary length, it will hit the light it was cast at. Therefore, you should
                // modify the time_bounds of your shadow ray to account for this. Using EPS_F is
                // recommended.

                Ray shadow_ray = Ray(hit.position, sample.direction, ray.lambda);
                shadow_ray.dist_bounds = Vec2(EPS_F, sample.distance);
                if(scene.hit(shadow_ray).hit) {
                    continue;
                }

                // Note: that along with the typical cos_theta, pdf factors, we divide by samples.
                // This is because we're doing another monte-carlo estimate of the lighting from
                // area lights here.
                radiance_out += (cos_theta / (samples * sample.pdf)) * Old2NewSpectrum(sample.radiance).sampleAtLambda(ray.lambda) * Old2NewSpectrum(attenuation).sampleAtLambda(ray.lambda);
            }
        };

        // If the BSDF is discrete (i.e. uses dirac deltas/if statements), then we are never
        // going to hit the exact right direction by sampling lights, so ignore them.
        if(!bsdf.is_discrete()) {

            // loop over all the lights and accumulate radiance.
            for(const auto& light : lights)
                sample_light(light);
            if(env_light.has_value())
                sample_light(env_light.value());
        }
    }

    // TODO (PathTracer): Task 5
    // Compute an indirect lighting estimate using path tracing with Monte Carlo.
    // (1) Ray objects have a depth field; if it reaches max_depth, you should
    // terminate the path.
    if(ray.depth >= max_depth) {
        return radiance_out;
    }

    // (2) Randomly select a new ray direction (it may be reflection or transmittance
    // ray depending on surface type) using bsdf.sample()

    BSDF_Sample bsdf_sample = bsdf.sample(out_dir, ray.lambda);

    // (3) Compute the throughput of the recursive ray. This should be the current ray's
    // throughput scaled by the BSDF attenuation, cos(theta), and BSDF sample PDF.
    // Potentially terminate the path using Russian roulette as a function of the new throughput.
    // Note that allowing the termination probability to approach 1 may cause extra speckling.

    float lambda_pdf = 1 /1.0f;
    Vec3 newDir = bsdf_sample.direction;
    float beta = bsdf_sample.attenuation * abs(bsdf_sample.direction.y) * (1 / (bsdf_sample.pdf*lambda_pdf));
    float recursive_ray_throughtput = beta * ray.throughput;
    // follow the sudo code on slide 58 from the "Global illumination" class

    float q = 1 - fmax(0, fmin(beta, 1));
    q *= 0.5;
    if(RNG::unit() < q){
        return radiance_out;
    }
    else {
        beta = beta / (1-q);
    }

    // (4) Create new scene-space ray and cast it to get incoming light. As with shadow rays,
    // you should modify time_bounds so that the ray does not intersect at time = 0. Remember to
    // set the new throughput and depth values.
    newDir = object_to_world.rotate(newDir);
    Ray scene_space_ray = Ray(hit.position, newDir,ray.lambda);
    scene_space_ray.throughput = recursive_ray_throughtput;
    scene_space_ray.depth = ray.depth + 1;
    scene_space_ray.dist_bounds = Vec2(EPS_F, std::numeric_limits<float>::infinity());

    float result = trace_ray(scene_space_ray);

    // (5) Add contribution due to incoming light with proper weighting. Remember to add in
    // the BSDF sample emissive term.
    radiance_out = radiance_out + result * beta + Old2NewSpectrum(bsdf_sample.emissive).sampleAtLambda(ray.lambda);
    return radiance_out;
}

} // namespace PT
