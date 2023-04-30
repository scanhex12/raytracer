#pragma once

#include <triangle.h>
#include <material.h>
#include <sphere.h>
#include <memory>

struct Object {
    std::shared_ptr<Material> material = nullptr;
    Triangle polygon;
    std::array<std::shared_ptr<Vector>, 3> normals;

    Object(std::shared_ptr<Material> mat, Vector a, Vector b, Vector c,
           std::shared_ptr<Vector> norm_a, std::shared_ptr<Vector> norm_b,
           std::shared_ptr<Vector> norm_c)
        : material(mat), polygon({a, b, c}), normals({norm_a, norm_b, norm_c}) {
    }

    const std::shared_ptr<Vector> GetNormal(size_t index) const {
        return normals[index];
    }
};

struct SphereObject {
    std::shared_ptr<Material> material = nullptr;
    Sphere sphere;

    SphereObject(std::shared_ptr<Material> mat, double x, double y, double z, double r)
        : material(mat), sphere(Vector({x, y, z}), r) {
    }
};
