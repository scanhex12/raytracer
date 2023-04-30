#pragma once

#include <vector.h>

class Triangle {
public:
    Triangle(std::initializer_list<Vector> list) {
        auto it = list.begin();
        for (size_t i = 0; i < 3; ++i) {
            vertices_[i] = *it;
            it++;
        }
    }
    double Area() const {
        Vector a = Difference(vertices_[0], vertices_[1]);
        Vector b = Difference(vertices_[0], vertices_[2]);
        return Length(CrossProduct(a, b)) / 2.0;
    }

    const Vector& GetVertex(size_t ind) const {
        return vertices_[ind];
    }

private:
    std::array<Vector, 3> vertices_;
};
