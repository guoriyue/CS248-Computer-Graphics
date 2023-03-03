
#include "../rays/bsdf.h"
#include "../util/rand.h"
#include "debug.h"

namespace PT {

Vec3 reflect(Vec3 dir) {

    // TODO (PathTracer): Task 6
    // Return reflection of dir about the surface normal (0,1,0).

    Vec3 surface_normal = Vec3(0.0f, 1.0f, 0.0f);
    Vec3 reflect_dir = -dir + 2.0f * dot(dir, surface_normal) * surface_normal;
    return reflect_dir;
}

Vec3 refract(Vec3 out_dir, float index_of_refraction, bool& was_internal) {

    // TODO (PathTracer): Task 6
    // Use Snell's Law to refract out_dir through the surface
    // Return the refracted direction. Set was_internal to false if
    // refraction does not occur due to total internal reflection,
    // and true otherwise.

    // When dot(out_dir,normal=(0,1,0)) is positive, then out_dir corresponds to a
    // ray exiting the surface into vaccum (ior = 1). However, note that
    // you should actually treat this case as _entering_ the surface, because
    // you want to compute the 'input' direction that would cause this output,
    // and to do so you can simply find the direction that out_dir would refract
    // _to_, as refraction is symmetric.
    Vec3 surface_normal = Vec3(0.0f, 1.0f, 0.0f);

    //check if the ray is same side as the normal
    float cos_theta_i = dot(surface_normal, out_dir);
    float eta_i_over_t;
    if(cos_theta_i >= 0) {
        eta_i_over_t = 1.0f/index_of_refraction;
    } else {
        eta_i_over_t = index_of_refraction/1.0f;
    }

    //check if the refraction is the total internal reflection
    float cos_theta_t_sq = 1 - pow(eta_i_over_t,2)*(1 - pow(cos_theta_i,2));
    if(cos_theta_t_sq < 0){
        was_internal = true;
        return Vec3();
    }

    was_internal = false;
    float cos_theta_t = sqrt(cos_theta_t_sq);
    if (cos_theta_i >= 0){
        cos_theta_t = -cos_theta_t;
    }
    
    Vec3 refract_dir = eta_i_over_t * (-out_dir) + (eta_i_over_t*dot(surface_normal, out_dir) + cos_theta_t )*surface_normal;
    return refract_dir;
}

BSDF_Sample BSDF_Lambertian::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 5
    // Implement lambertian BSDF. Use of BSDF_Lambertian::sampler may be useful

    BSDF_Sample ret;
    Vec3 sampled_direction = sampler.sample(ret.pdf); // What was the PDF of the sampled direction?
    ret.attenuation = evaluate(
        out_dir, sampled_direction); // What is the ratio of reflected/incoming light? The higher the brighter
    ret.direction = sampled_direction; // What direction should we sample incoming light from?
    return ret;
}

Spectrum BSDF_Lambertian::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    return albedo * (1.0f / PI_F);
}

BSDF_Sample BSDF_Mirror::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    // Implement mirror BSDF

    BSDF_Sample ret;
    Vec3 surface_normal = Vec3(0.0f, 1.0f, 0.0f);
    ret.attenuation = reflectance/ dot(out_dir, surface_normal);          // What is the ratio of reflected/incoming light?
    ret.direction = reflect(out_dir);       // What direction should we sample incoming light from?
    ret.pdf = 1.0f; // Was was the PDF of the sampled direction? (In this case, the PMF)
    return ret;
}

Spectrum BSDF_Mirror::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // Technically, we would return the proper reflectance
    // if in_dir was the perfectly reflected out_dir, but given
    // that we assume these are single exact directions in a
    // continuous space, just assume that we never hit them
    // _exactly_ and always return 0.
    return {};
}

BSDF_Sample BSDF_Glass::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6

    // Implement glass BSDF.
    // (1) Compute Fresnel coefficient. Tip: use Schlick's approximation.
    // (2) Reflect or refract probabilistically based on Fresnel coefficient. Tip: RNG::coin_flip
    // (3) Compute attenuation based on reflectance or transmittance

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?

    BSDF_Sample ret;
    ret.attenuation = Spectrum(); // What is the ratio of reflected/incoming light?
    ret.direction = Vec3();       // What direction should we sample incoming light from?
    ret.pdf = 0.0f; // Was was the PDF of the sampled direction? (In this case, the PMF)

    //Schlick's approximation
    float R0 = pow(((1.0f - index_of_refraction)/(1+ index_of_refraction)),2);
    Vec3 surface_normal = Vec3(0.0f, 1.0f, 0.0f);
    float cos_theta = abs(dot(out_dir, surface_normal));
    float R_theta = R0 + (1-R0) * pow((1-cos_theta),5);

    // check if it is total internal reflection
    bool was_internal;
    Vec3 refract_dir = refract(out_dir, index_of_refraction, was_internal);

    if (was_internal){
        ret.attenuation = R_theta * reflectance / cos_theta;
        ret.direction = reflect(out_dir);
        ret.pdf = R_theta;
        return ret;
    }
    else {
        //not total internal reflection
        // (e.g., If the Fresnel reflectance is 0.9, then you should generate a reflection ray 90% of the time
        bool is_reflection = RNG::coin_flip(R_theta);
        if(is_reflection) {
            ret.attenuation = R_theta * reflectance / cos_theta;
            ret.direction = reflect(out_dir);
            ret.pdf = R_theta;
            return ret;
        }
        ret.attenuation = (1.0f - R_theta) * transmittance / cos_theta;
        ret.direction = refract_dir;
        ret.pdf = 1.0f - R_theta;
        return ret;
    }

    
}

Spectrum BSDF_Glass::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // As with BSDF_Mirror, just assume that we never hit the correct
    // directions _exactly_ and always return 0.
    return {};
}

BSDF_Sample BSDF_Diffuse::sample(Vec3 out_dir) const {
    BSDF_Sample ret;
    ret.direction = sampler.sample(ret.pdf);
    ret.emissive = radiance;
    ret.attenuation = {};
    return ret;
}

Spectrum BSDF_Diffuse::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // No incoming light is reflected; only emitted
    return {};
}

BSDF_Sample BSDF_Refract::sample(Vec3 out_dir) const {

    // TODO (PathTracer): Task 6
    // Implement pure refraction BSDF.

    // Be wary of your eta1/eta2 ratio - are you entering or leaving the surface?

    BSDF_Sample ret;

    bool was_internal;
    Vec3 surface_normal = Vec3(0.0f, 1.0f, 0.0f);
    ret.attenuation = transmittance / abs(dot(out_dir, surface_normal)); // What is the ratio of reflected/incoming light?
    ret.direction = refract(out_dir, index_of_refraction, was_internal);       // What direction should we sample incoming light from?
    ret.pdf = 1.0f; // Was was the PDF of the sampled direction? (In this case, the PMF)
    return ret;
}

Spectrum BSDF_Refract::evaluate(Vec3 out_dir, Vec3 in_dir) const {
    // As with BSDF_Mirror, just assume that we never hit the correct
    // directions _exactly_ and always return 0.
    return {};
}

} // namespace PT
