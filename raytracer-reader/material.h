#pragma once

#include <vector.h>
#include <string>

struct Material {
    std::string name;
    Vector ambient_color;
    Vector diffuse_color;
    Vector specular_color;
    Vector intensity;
    double specular_exponent;
    double refraction_index;
    std::array<double, 3> albedo;
    double tr;

    Material() {
        albedo[0] = 1;
    }

    Material(std::string mname, Vector ka, Vector kd, Vector ks, Vector ke, double ns, double ni,
             std::array<double, 3> al, double tr)
        : name(mname),
          ambient_color(ka),
          diffuse_color(kd),
          specular_color(ks),
          intensity(ke),
          specular_exponent(ns),
          refraction_index(ni),
          tr(tr) {
        albedo[0] = al[0];
        albedo[1] = al[1];
        albedo[2] = al[2];
    }
};
