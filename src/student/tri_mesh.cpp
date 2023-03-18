#include "../rays/tri_mesh.h"
#include "debug.h"

namespace PT {

BBox Triangle::bbox() const {

    // TODO (PathTracer): Task 2
    // compute the bounding box of the triangle

    // Beware of flat/zero-volume boxes! You may need to
    // account for that here, or later on in BBox::intersect

    BBox box;

    float x_min = std::min(vertex_list[v0].position.x, std::min(vertex_list[v1].position.x, vertex_list[v2].position.x));
    float x_max = std::max(vertex_list[v0].position.x, std::max(vertex_list[v1].position.x, vertex_list[v2].position.x));
    float y_min = std::min(vertex_list[v0].position.y, std::min(vertex_list[v1].position.y, vertex_list[v2].position.y));
    float y_max = std::max(vertex_list[v0].position.y, std::max(vertex_list[v1].position.y, vertex_list[v2].position.y));
    float z_min = std::min(vertex_list[v0].position.z, std::min(vertex_list[v1].position.z, vertex_list[v2].position.z));
    float z_max = std::max(vertex_list[v0].position.z, std::max(vertex_list[v1].position.z, vertex_list[v2].position.z));
    
    box.min.x = x_min;
    box.max.x = x_max;
    box.min.y = y_min;
    box.max.y = y_max;
    box.min.z = z_min;
    box.max.z = z_max;

    return box;
}

Trace Triangle::hit(const Ray& ray) const {
    // printf("Triangle::hit\n");
    // Vertices of triangle - has postion and surface normal
    // See rays/tri_mesh.h for a description of this struct
    
    Tri_Mesh_Vert v_0 = vertex_list[v0];
    Tri_Mesh_Vert v_1 = vertex_list[v1];
    Tri_Mesh_Vert v_2 = vertex_list[v2];

    // here just to avoid unused variable warnings, students should remove the following three lines.
    (void)v_0;
    (void)v_1;
    (void)v_2;
    
    // TODO (PathTracer): Task 2
    // Intersect this ray with a triangle defined by the above three points.
    // Intersection should yield a ray t-value, and a hit point (u,v) on the surface of the triangle

    // You'll need to fill in a "Trace" struct describing information about the hit (or lack of hit)

    Trace ret;
    ret.origin = ray.point;
    ret.hit = false;       // was there an intersection?
    ret.distance = 0.0f;   // at what distance did the intersection occur?
    ret.position = Vec3{}; // where was the intersection?
    ret.normal = Vec3{};   // what was the surface normal at the intersection?
                           // (this should be interpolated between the three vertex normals)

    Vec3 e_1 = v_1.position - v_0.position;
    Vec3 e_2 = v_2.position - v_0.position;
    Vec3 s = ray.point - v_0.position;
    auto determinant = [](auto e_1, auto e_2, auto dir) {
        return dot(cross(e_1, dir), e_2);
    };

    // float denominator = dot(cross(e_1, ray.dir), e_2);
    float denominator = determinant(e_1, e_2, ray.dir);
    if (std::abs(denominator) < 1e-6) {
        // parallel to the plane
        return ret;
    }
    Vec3 uvt = Vec3(-1 * determinant(s, ray.dir, e_2), determinant(e_1, s, ray.dir), -1 * determinant(s, e_1, e_2)) / denominator;
    
    float u = uvt[0];
    float v = uvt[1];
    float t = uvt[2];
    

    if (0 <= u && u <= 1 && 0 <= v && v <= 1 && u + v <= 1 && ray.dist_bounds[0] <= t && t <= ray.dist_bounds[1]) {
        // printf("u: %f, v: %f, t: %f\n", u, v, t);
        // printf("ray.dist_bounds[0] %f, ray.dist_bounds[1] %f\n", ray.dist_bounds[0], ray.dist_bounds[1]);
        ret.hit = true;
        ret.distance = t;
        ret.position = ray.point + ray.dir * t;
        ret.normal = (1 - u - v) * v_0.normal + u * v_1.normal + v * v_2.normal;
    }

    return ret;
}

Triangle::Triangle(Tri_Mesh_Vert* verts, unsigned int v0, unsigned int v1, unsigned int v2)
    : vertex_list(verts), v0(v0), v1(v1), v2(v2) {
}

void Tri_Mesh::build(const GL::Mesh& mesh) {

    verts.clear();
    triangles.clear();

    for(const auto& v : mesh.verts()) {
        verts.push_back({v.pos, v.norm});
    }

    const auto& idxs = mesh.indices();

    std::vector<Triangle> tris;
    for(size_t i = 0; i < idxs.size(); i += 3) {
        tris.push_back(Triangle(verts.data(), idxs[i], idxs[i + 1], idxs[i + 2]));
    }

    triangles.build(std::move(tris), 4);
}

Tri_Mesh::Tri_Mesh(const GL::Mesh& mesh) {
    build(mesh);
}

Tri_Mesh Tri_Mesh::copy() const {
    Tri_Mesh ret;
    ret.verts = verts;
    ret.triangles = triangles.copy();
    return ret;
}

BBox Tri_Mesh::bbox() const {
    return triangles.bbox();
}

Trace Tri_Mesh::hit(const Ray& ray) const {
    Trace t = triangles.hit(ray);
    return t;
}

size_t Tri_Mesh::visualize(GL::Lines& lines, GL::Lines& active, size_t level,
                           const Mat4& trans) const {
    return triangles.visualize(lines, active, level, trans);
}

} // namespace PT
