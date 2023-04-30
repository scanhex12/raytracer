#pragma once

#include <vector.h>

class Ray {
public:
    Ray(Vector origin, Vector direction) : origin_(origin), direction_(direction) {
    }
    const Vector& GetOrigin() const {
        return origin_;
    }
    const Vector& GetDirection() const {
        return direction_;
    }

    void Step(double eps) {
        Vector norm = direction_;
        norm.Normalize();
        origin_ = Sum(origin_, Mult(norm, eps));
    }
private:
    Vector origin_;
    Vector direction_;
};
