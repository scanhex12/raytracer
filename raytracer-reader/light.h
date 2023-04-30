#pragma once

#include <vector.h>
#include <triangle.h>

struct Light {
    Vector position;
    Vector intensity;

    Light(Vector position, Vector intensity)
        : position(std::move(position)), intensity(std::move(intensity)) {
    }

    Light(double x, double y, double z, double r, double g, double b)
        : position(Vector({x, y, z})), intensity(Vector({r, g, b})) {
    }
};
